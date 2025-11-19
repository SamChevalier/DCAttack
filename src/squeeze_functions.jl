function zero_nonlinear_costs!(basic_network_data)
    for (key,val) in basic_network_data["gen"]
        val["cost"][1] = 0.0
        val["cost"][3] = 0.0
    end
end

function solve_dcopf(basic_network_data)
    # pm_result = PowerModels.solve_dc_opf(basic_network_data, Ipopt.Optimizer)
    pm_result = PowerModels.solve_dc_opf(basic_network_data, Gurobi.Optimizer)

    # system information
    nb = length(basic_network_data["bus"])
    nd = length(basic_network_data["load"])
    nl = length(basic_network_data["branch"])

    # how many gens can produce power AND have non-zero flexibility?
    ng = 0
    for (key,val) in basic_network_data["gen"]
        if val["pmax"] > 0.0
            if !isapprox(val["pmax"],val["pmin"])
                ng += 1
            end
        end
    end

    # map buses to indices
    index = 1
    bus_id_list    = Int64[]
    bus_index_list = Int64[]
    for (key,val) in basic_network_data["bus"]
        # => don't use "source_id"
        source_id = val["bus_i"]
        push!(bus_id_list, source_id)
        push!(bus_index_list, index)
        index += 1
    end

    function map_id_to_index(bus_id, bus_id_list)
        findall(bus_id_list .== bus_id)[1]
    end

    # get the loading at each bus!
    pd = zeros(nb)
    pert_index = []
    for (key,val) in basic_network_data["load"]
        # don't use "source_id"
        bus_id = val["load_bus"]
        index = map_id_to_index(bus_id, bus_id_list)
        pd[index] += val["pd"]

        # is there nonzero active power loading at this bus?
        if abs(val["pd"]) > 0.0
            push!(pert_index, index)
        end
    end

    # get the generation mapping: N (nb by ng)
    N        = zeros(nb,ng)
    pmax     = zeros(ng)
    pmin     = zeros(ng)
    lin_cost = zeros(ng)
    Cost1    = zeros(ng)
    Cost2    = zeros(ng)
    Cost3    = zeros(ng)
    Pg_soln  = zeros(ng)

    ii = 1
    slack_idx = 0
    slack_gen = 0
    c0        = 0.0 # this is the cost of non-flexible generators
    for (key,val) in basic_network_data["gen"]
        if val["pmax"] > 0.0
            if !isapprox(val["pmax"],val["pmin"])
                # don't use "source_id"
                gen_id = val["gen_bus"]
                index = map_id_to_index(gen_id, bus_id_list)

                pmax[ii]  = val["pmax"]
                pmin[ii]  = val["pmin"]
                Cost1[ii] = val["cost"][1]
                Cost2[ii] = val["cost"][2]
                Cost3[ii] = val["cost"][3]
                N[index, ii] = 1
                Pg_soln[ii]  = pm_result["solution"]["gen"][key]["pg"]

                # also, choose the slack generator here
                if slack_idx == 0
                    slack_idx = index
                    slack_gen = ii
                end

                ii += 1
            else
                # in this case, just add the generation to the load
                gen_id = val["gen_bus"]
                index = map_id_to_index(gen_id, bus_id_list)
                pd[index] -= val["pmax"]
                c0 += val["pmax"]*val["cost"][2]
            end
        end
    end
    lin_cost = Cost2

    # branch information
    b_line_list = Float64[]
    to_bus_list = Int64[]
    fr_bus_list = Int64[]
    flow_max    = Float64[]
    E = zeros(nl,nb)
    ii = 1

    for (key,val) in basic_network_data["branch"]
        fbus = map_id_to_index(val["f_bus"], bus_id_list)
        tbus = map_id_to_index(val["t_bus"], bus_id_list)
        push!(fr_bus_list, fbus)
        push!(to_bus_list, tbus)
        push!(b_line_list, imag(1/(val["br_r"] + im*val["br_x"])))
        push!(flow_max, val["rate_c"])
        E[ii,fbus] = -1
        E[ii,tbus] = 1
        ii = ii + 1
    end

    # build ptdf
    Ehat = E[:,Not(slack_idx)]
    Yl   = diagm(b_line_list)
    ptdf = Yl*Ehat*inv(Ehat'*Yl*Ehat)
    ptdf_full = [ptdf[:,1:(slack_idx-1)] zeros(nl) ptdf[:,slack_idx:end]]

    # lightly sparsify...
    ptdf_full[ abs.(ptdf_full) .< 1e-6 ] .= 0.0

    # Build model
    model = Model(Gurobi.Optimizer)
    @variable(model, Pg[1:ng])

    # Extract A, B, C, D
    A_lp = ones(Float64, ng)'
    univec = ones(length(pd))
    B_lp = -[univec'*pd]
    C_lp =  [ptdf_full*N;-ptdf_full*N;I(ng);-I(ng)]
    D_lp = -[ptdf_full*pd.+flow_max;-ptdf_full*pd.+flow_max;pmax;-pmin]

    # solve with Gurobi using PTDF
    c1 = @constraint(model,  A_lp*Pg  .+ B_lp[1] .== 0.0)
    c2 = @constraint(model , C_lp*Pg  .+ D_lp    .<= 0.0)
    # => quad objective: o = sum(Cost1[i]*(Mσ[i,i]*Pg[i] + mμ[i]).^2 + Cost2[i]*(Mσ[i,i]*Pg[i] + mμ[i]) + Cost3[i] for i in 1:length(Pg))
    @objective(model, Min, lin_cost'*Pg + c0)
    @time optimize!(model)
    println(objective_value(model))
    objective_value(model)

    # as a final step, build the inequality-based linear system
    ptdf_reduced = copy(ptdf)
    println(slack_idx)
    Ng           = N[Not(slack_idx), Not(slack_gen)]
    n_perts      = length(pert_index)
    N_delta      = zeros(nb,n_perts)
    for ii in 1:n_perts
        N_delta[pert_index[ii],ii] = 1
    end

    g_sum = ones(ng-1)
    d_sum = ones(n_perts)

    # define the paper matrices
    A = [ptdf_reduced*Ng; 
        -ptdf_reduced*Ng; 
        -g_sum'; 
        Matrix(I,ng-1,ng-1); 
        g_sum'; 
        -Matrix(I,ng-1,ng-1)]
    B = [-ptdf_full*N_delta; ptdf_full*N_delta; d_sum'; zeros(ng-1,n_perts); -d_sum'; zeros(ng-1,n_perts)]
    b = [-ptdf_full*pd - flow_max; -flow_max + ptdf_full*pd; [sum(pd); zeros(ng-1)] - pmax; pmin - [sum(pd); zeros(ng-1)]]

    # as a final test, make sure this model is also correct
    model_fl = Model(Gurobi.Optimizer)
    @variable(model_fl, pg[2:ng])
    @constraint(model_fl, A*pg + b .<= 0.0)
    @objective(model_fl, Min, lin_cost[2:end]'*pg + lin_cost[1]*(sum(pd)-sum(pg)) + c0)
    optimize!(model_fl)

    # pack it all into a dictionary
    sys = Dict(:pd => pd,
               :lin_cost => lin_cost,
               :univec => univec, 
               :flow_max => flow_max, 
               :pmax => pmax, 
               :pmin => pmin, 
               :N => N, 
               :N_delta => N_delta,
               :ptdf_full => ptdf_full, 
               :ptdf => ptdf, 
               :A_lp => A_lp, 
               :B_lp => B_lp, 
               :C_lp => C_lp, 
               :D_lp => D_lp, 
               :nb => nb, 
               :ng => ng, 
               :nl => nl, 
               :pert_index => pert_index, 
               :A => A, 
               :B => B, 
               :b => b,
               :c0 => c0)

    # output
    return pm_result, model, model_fl, sys
end

function solve_dcopf_ABCD(A, B, C, D, pmax, pmin, lin_cost)
    model = Model(Gurobi.Optimizer)
    @variable(model, pg[1:ng])
    @constraint(model, A*pg .+ B .== 0.0)
    @constraint(model, C*pg .+ D .<= 0.0)
    @objective(model, Min, lin_cost'*pg)

    @time optimize!(model)
    println("Optimal g: ", JuMP.value.(pg[1:ng]))
    println("Total lin_cost: ", objective_value(model))

    return model
end

function test_parsing(cases::Vector{String})
    case_status1 = []
    case_status2 = []
    for case in cases
        network_data       = pglib(case)
        basic_network_data = PowerModels.make_basic_network(network_data)
        zero_nonlinear_costs!(basic_network_data)

        # solve dcopf
        pm_result, model, model_fl, sys = solve_dcopf(basic_network_data);
        println()
        println(pm_result["objective"])
        println(objective_value(model))
        push!(case_status1, isapprox(pm_result["objective"],objective_value(model)))
        push!(case_status2, isapprox(pm_result["objective"],objective_value(model_fl)))
    end
    println()
    println("==========================================")
    println()
    
    if all([case_status1; case_status2])
        println("All cases passed ✔")
    else
        println("At least one case failed ✖")
        println([case_status1 case_status2])
    end
end

function solve_inequality_dcopf(A, B, C, D, lin_cost)
    model = Model(Gurobi.Optimizer)
    @variable(model, pg[2:ng])
    @constraint(model, (C[:,2:end]-C[:,1]*A[2:end]')*pg .+ (D + (C[:,1]*.-B[1])) .<= 0.0)
    @objective(model, Min, lin_cost[2:end]'*pg + lin_cost[1]*(-A[2:end]'*pg-B[1]))

    @time optimize!(model)
    println(objective_value(model))

    # build structs
    C_in = C[:,2:end]-C[:,1]*A[2:end]'
    d_in = D + (C[:,1]*.-B[1])

    return model, C_in, d_in
end

# feasibility model
function dcfeas_test(pd,A,univec,ptdf_full,N,ng,flow_max,pmax,pmin)
    B = -[univec'*pd]
    C =  [ptdf_full*N;-ptdf_full*N;I(ng);-I(ng)]
    D = -[ptdf_full*pd.+flow_max;-ptdf_full*pd.+flow_max;pmax;-pmin]

    model = Model(Gurobi.Optimizer)
    @variable(model, Pg[1:ng])
    @constraint(model, A*Pg .+ B .== 0.0)
    @constraint(model, C*Pg .+ D .<= 0.0)
    @objective(model, Min, 0.0)

    optimize!(model)

    return Int(termination_status(model)), model
end

# feasibility model -- inequality!
function dcfeas_ineq(C_in, d_in, gamma)

    model = Model(Gurobi.Optimizer)
    @variable(model, x[1:size(C_in,2)])
    @variable(model, s)
    @constraint(model, 0.0 .<= s)
    @constraint(model, C_in*x .+ d_in .+ s .<= -gamma)
    @objective(model, Max, s)

    optimize!(model)

    return Int(termination_status(model)), value.(x)
end

function solve_farkas_lemma(A, B, b, MIPGap, TimeLimit, nb, x0; init=true, extra_string="")
    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "MIPGap", MIPGap)
    set_optimizer_attribute(model, "TimeLimit", TimeLimit)
    set_optimizer_attribute(model, "MIPFocus", 3)

    #set_optimizer_attribute(model, "OptimalityTol", 1e-5)
    #set_optimizer_attribute(model, "FeasibilityTol", 1e-5)

    num_mu  = size(A,1)
    n_perts = size(B,2)
    @variable(model, delta[1:n_perts])
    @variable(model, mu[1:num_mu])
    @constraint(model, 0.0 .<= mu)
    @constraint(model, A'*mu .== 0.0)
    @constraint(model, mu'*(B*delta+b) .== 0.001)
    @objective(model, Min, dot(delta,delta))

    # initialize with x0
    if init
        set_start_value.(delta, x0[:delta])
        set_start_value.(mu, x0[:mu])
    end

    # log structures
    global timelog  = []
    global bestlog  = []
    global boundlog = []
    MOI.set(model, Gurobi.CallbackFunction(), callback_log)

    # optimize
    optimize!(model)

    # clean up
    faraks_log = filterlog(timelog, bestlog, boundlog)

    # write to file
    testdata_file = "./data/farkas_"*string(nb)*"bus"*extra_string*".h5"
    write_hdf5data(testdata_file, faraks_log)

    return faraks_log
end

function solve_farkas_lemma_record_log(A, B, b, MIPGap, TimeLimit, nb, x0; init=true, extra_string="")
    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "MIPGap", MIPGap)
    set_optimizer_attribute(model, "TimeLimit", TimeLimit)
    set_optimizer_attribute(model, "MIPFocus", 3)
    set_optimizer_attribute(model, "LogFile", "./data/gurobi_branching_"*string(nb)*"bus"*".log")

    num_mu  = size(A,1)
    n_perts = size(B,2)
    @variable(model, delta[1:n_perts])
    @variable(model, mu[1:num_mu])
    @constraint(model, 0.0 .<= mu)
    @constraint(model, A'*mu .== 0.0)
    @constraint(model, mu'*(B*delta+b) .== 0.001)
    @objective(model, Min, dot(delta,delta))

    # initialize with x0
    if init
        set_start_value.(delta, x0[:delta])
        set_start_value.(mu, x0[:mu])
    end

    # log structures
    global timelog  = []
    global bestlog  = []
    global boundlog = []
    MOI.set(model, Gurobi.CallbackFunction(), callback_log)

    # optimize
    optimize!(model)

    # clean up
    # => faraks_log = filterlog(timelog, bestlog, boundlog)
    # write to file
    # => testdata_file = "./farkas_"*string(nb)*"bus"*extra_string*".h5"
    # => write_hdf5data(testdata_file, faraks_log)
    # =>return faraks_log
end

function solve_farkas_lemma_local(A, B, b; initialize_start = false, lower_ipopt_tol = false, delta0, mu0)
    model = Model(Ipopt.Optimizer)

    if lower_ipopt_tol == true
        set_optimizer_attribute(model, "tol", 1e-5)
        set_optimizer_attribute(model, "constr_viol_tol", 1e-5)
        set_optimizer_attribute(model, "dual_inf_tol", 1e-5)
    end

    num_mu  = size(A,1)
    n_perts = size(B,2)
    @variable(model, delta[1:n_perts])
    @variable(model, mu[1:num_mu])
    @constraint(model, 0.0 .<= mu)
    @constraint(model, A'*mu .== 0.0)
    @constraint(model, mu'*(B*delta+b) .== 0.001)
    @objective(model, Min, dot(delta,delta))

    if initialize_start == true
        set_start_value.(delta, delta0)
        set_start_value.(mu,    mu0)
    end

    # optimize
    optimize!(model)

    # status?
    println(termination_status(model))

    return termination_status(model), value.(delta), value.(mu), value(dot(delta,delta))
end

function solve_control_local(A, B, b, pg0, t0, v0, tmax; random_start=false)
    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "max_cpu_time", tmax)

    set_optimizer_attribute(model, "tol", 1e-6)
    set_optimizer_attribute(model, "constr_viol_tol", 1e-6)
    set_optimizer_attribute(model, "dual_inf_tol", 1e-6)

    # add the control policy
    nc = size(A,1)
    nx = size(A,2)
    n_perts = size(B,2)
    @variable(model, t)
    @variable(model, G[1:nx,1:n_perts])
    @variable(model, p0[1:nx])
    @variable(model, v[1:nc])

    if random_start == true
        set_start_value.(G, 0.0 .+ randn(nx, n_perts))
        set_start_value.(p0, pg0 + 0.1*randn(nx))
        set_start_value(t, t0 + 0.1*randn())
        set_start_value.(v, v0 + 0.1*randn(nc))
    else
        set_start_value.(G, 0.0)
        set_start_value.(p0, pg0)
        set_start_value(t, t0)
        set_start_value.(v, v0)
    end

    # constrain
    for ii in 1:nc
        @constraint(model, v[ii] == ((G'*A[ii,:] + B[ii,:])'*(G'*A[ii,:] + B[ii,:])))
        @constraint(model, t*v[ii] .<= (A[ii,:]'*p0 + b[ii])^2)
    end

    # base feasibility constraint
    @constraint(model, A*p0 + b .<= 0.0)
    @objective(model, Max, t)

    # optimize
    optimize!(model)

    # status?
    println(termination_status(model))

    return termination_status(model), value.(G), value.(p0), value.(v), value(t)
end

function solve_control(A, B, b, MIPGap, TimeLimit, nb, x0; init=true, extra_string="")
    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "MIPGap", MIPGap)
    set_optimizer_attribute(model, "TimeLimit", TimeLimit)
    set_optimizer_attribute(model, "MIPFocus", 3)

    set_optimizer_attribute(model, "OptimalityTol", 1e-5)
    set_optimizer_attribute(model, "FeasibilityTol", 1e-5)

    # add the control policy
    nc = size(A,1)
    nx = size(A,2)
    n_perts = size(B,2)
    @variable(model, t)
    @variable(model, G[1:nx,1:n_perts])
    @variable(model, p0[1:nx])
    @variable(model, v[1:nc])

    if init
        set_start_value.(G, x0[:G])
        set_start_value.(p0, x0[:pg0])
        set_start_value.(v,x0[:v])
        set_start_value(t, x0[:t])
    end

    # constrain
    for ii in 1:nc
        @constraint(model, v[ii] == ((G'*A[ii,:] + B[ii,:])'*(G'*A[ii,:] + B[ii,:])))
        @constraint(model, t*v[ii] .<= (A[ii,:]'*p0 + b[ii])^2)
    end

    # base feasibility constraint
    @constraint(model, A*p0 + b .<= 0.0)
    @objective(model, Max, t)

    # log structures
    global timelog  = []
    global bestlog  = []
    global boundlog = []
    MOI.set(model, Gurobi.CallbackFunction(), callback_log)

    # optimize
    optimize!(model)

    # status?
    control_log = filterlog(timelog, bestlog, boundlog)

    # write to file
    testdata_file = "./data/control_"*string(nb)*"bus"*extra_string*".h5"
    write_hdf5data(testdata_file, control_log)

    return control_log
end



function solve_squeeze(A, B, b, MIPGap, TimeLimit, nb, x0)
    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "MIPGap", MIPGap)
    set_optimizer_attribute(model, "TimeLimit", TimeLimit)

    num_mu  = size(A,1)
    n_perts = size(B,2)
    @variable(model, delta[1:n_perts])
    @variable(model, mu[1:num_mu])
    @constraint(model, 0.0 .<= mu)
    @constraint(model, A'*mu .== 0.0)
    @constraint(model, mu'*(B*delta+b) .== 0.001)

    # initialize with x0
    set_start_value.(delta, x0[:delta])
    set_start_value.(mu, x0[:mu])

    # add the control policy
    nc = size(A,1)
    nx = size(A,2)
    @variable(model, t)
    @variable(model, G[1:nx,1:n_perts])
    @variable(model, p0[1:nx])
    @variable(model, v[1:nc])

    # initialize
    set_start_value.(G,  x0[:G])
    set_start_value.(p0, x0[:pg0])
    set_start_value.(v,  x0[:v])
    set_start_value.(t,  x0[:t])


    for ii in 1:nc
        @constraint(model, v[ii] == ((G'*A[ii,:] + B[ii,:])'*(G'*A[ii,:] + B[ii,:])))
        @constraint(model, t*v[ii] .<= (A[ii,:]'*p0 + b[ii])^2)
    end

    # base feasibility constraint
    @constraint(model, A*p0 + b .<= 0.0)
    @objective(model, Min,  (dot(delta,delta) - t)^2)

    # log structures
    global timelog  = []
    global bestlog  = []
    global boundlog = []
    MOI.set(model, Gurobi.CallbackFunction(), callback_log)

    # optimize
    optimize!(model)

    # clean up
    squeeze_log = filterlog(timelog, bestlog, boundlog)

    # write to file
    testdata_file = "./data/squeeze_"*string(nb)*"bus.h5"
    t_final       = value(t)
    delta2_final  = value(dot(delta,delta))
    write_squeeze_hdf5data(testdata_file, squeeze_log, t_final, delta2_final)

    return squeeze_log
end


# Let's try to build Farkas' Lemma now
function farkas_test(pd,A,univec,ptdf_full,N,ng,flow_max,pmax,pmin)
    @warn("depricated")
    B = -univec'*pd
    C =  [ptdf_full*N;-ptdf_full*N;I(ng);-I(ng)]
    D = -[ptdf_full*pd.+flow_max;-ptdf_full*pd.+flow_max;pmax;-pmin]

    model = Model(Gurobi.Optimizer)
    @variable(model, lambda)
    nm = size(C,1)
    @variable(model, 0.0 <= mu[1:nm])
    @constraint(model, lambda*A + mu'*C .== 0)
    @constraint(model, 1.0 .<= lambda*B + mu'*D)
    @objective(model, Min, 0)
    optimize!(model)

    return Int(termination_status(model))
end

function find_smallest_infeasibility_perturbation(nb, ng, A, C, univec, pd, ptdf_full, flow_max, pmax, pmin, zero_perts)
    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "MIPGap", 0.01)
    set_optimizer_attribute(model, "OutputFlag", 1)    # Enable output
    set_optimizer_attribute(model, "LogToConsole", 1)  # Ensure output is printed to terminal
    set_optimizer_attribute(model, "LogFile", "gurobi_branching.log")
    set_optimizer_attribute(model, "DisplayInterval", 1)   # Log every node

    @variable(model, gamma[1:nb])
    @constraint(model, gamma[1:zero_perts] .== 0.0)

    # @constraint(model, dot(gamma,gamma) .<= 0.57)

    B = -univec'*(pd .+ gamma[1:nb])
    D = -[ptdf_full*(pd .+ gamma[1:nb]).+flow_max;-ptdf_full*(pd .+ gamma).+flow_max;pmax;-pmin]

    @variable(model, lambda)
    nm = size(C,1)
    @variable(model, 0.0 <= mu[1:nm])
    @constraint(model, lambda*A + mu'*C .== 0)
    @constraint(model, 0.001 ==  lambda*B + mu'*D)
    # => @constraint(model, 0.001 .<=  lambda*B + mu'*D)
    @objective(model, Min, dot(gamma,gamma))
    optimize!(model)

    println("termination status: ", termination_status(model))
    return model
end

function find_smallest_infeasibility_perturbation_explicit_comparison(nb, ng, A, C, univec, pd, ptdf_full, flow_max, pmax, pmin, zero_perts)
    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "MIPGap", 0.01)
    @variable(model, gamma[1:nb])
    @constraint(model, gamma[1:zero_perts] .== 0.0)

    @variable(model, lambda)
    nm = size(C,1)
    @variable(model, 0.0 <= mu[1:nm])

    @constraint(model, lambda*A + mu'*C .== 0)

    mu1 = mu[(1):(nl)]
    mu2 = mu[(nl+1):(2*nl)]
    mu3 = mu[(2*nl+1):(2*nl+ng)]
    mu4 = mu[(2*nl+ng+1):(2*nl+2*ng)]

    t1 = -lambda*dot(univec,pd)
    t2 = -lambda*dot(univec,gamma)

    t3 = -mu1'*ptdf_full*pd
    t4 = -mu1'*ptdf_full*gamma
    t5 = -mu1'*flow_max

    t6 = mu2'*ptdf_full*pd
    t7 = mu2'*ptdf_full*gamma 
    t8 = -mu2'*flow_max

    t9 = -mu3'*pmax
    t10= +mu4'*pmin

    @constraint(model, 0.001 == t1+t2+t3+t4+t5+t6+t7+t8+t9+t10)
    @objective(model, Min, dot(gamma,gamma))
    optimize!(model)

    termination_status(model)
    println(value.(gamma))

    # build PSD test
    v = [1; value(lambda); value.(mu); value.(gamma)]
    M = v*(v')

    lam_idx = 2
    mu_idx  = (2+1):(2 + 2*nl+2*ng)
    gam_idx = (2 + 2*nl+2*ng + 1):(2 + 2*nl+2*ng + nb)

    mu1_idx = mu_idx[(1):(nl)]
    mu2_idx = mu_idx[(nl+1):(2*nl)]
    mu3_idx = mu_idx[(2*nl+1):(2*nl+ng)]
    mu4_idx = mu_idx[(2*nl+ng+1):(2*nl+2*ng)]

    lambda_n  = M[lam_idx,1]
    mu_n      = M[mu_idx,1]
    gamma_n   = M[gam_idx,1]

    mu1_n = mu_n[(1):(nl)]
    mu2_n = mu_n[(nl+1):(2*nl)]
    mu3_n = mu_n[(2*nl+1):(2*nl+ng)]
    mu4_n = mu_n[(2*nl+ng+1):(2*nl+2*ng)]

    lam_gam_n = M[lam_idx,gam_idx]
    gam_mu1_n = M[gam_idx,mu1_idx]
    gam_mu2_n = M[gam_idx,mu2_idx]

    # build the t's and compare
    t1n  = -lambda_n*dot(univec,pd)
    t2n  = -sum(lam_gam_n)
    t3n  = -mu1_n'*ptdf_full*pd
    t4n  = -tr(ptdf_full*gam_mu1_n)
    t5n  = -mu1_n'*flow_max
    t6n  = mu2_n'*ptdf_full*pd
    t7n  = tr(ptdf_full*gam_mu2_n) 
    t8n  = -mu2_n'*flow_max
    t9n  = -mu3_n'*pmax
    t10n = +mu4_n'*pmin

    return model
end

function find_smallest_infeasibility_perturbation_McCormick(nb, ng, A, C, univec, pd, ptdf_full, flow_max, pmax, pmin, zero_perts, gamma_min, gamma_max, lambda_min, lambda_max, mu1_min, mu1_max, mu2_min, mu2_max)
    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "MIPGap", 0.01)
    @variable(model, gamma[1:nb])
    @constraint(model, gamma[1:zero_perts] .== 0.0)

    @constraint(model, dot(gamma,gamma) .<= 0.57)

    @variable(model, lambda)
    nm = size(C,1)
    @variable(model, 0.0 <= mu[1:nm])
    @constraint(model, lambda*A + mu'*C .== 0)

    mu1 = mu[(1):(nl)]
    mu2 = mu[(nl+1):(2*nl)]
    mu3 = mu[(2*nl+1):(2*nl+ng)]
    mu4 = mu[(2*nl+ng+1):(2*nl+2*ng)]

    # add intermediate variables
    delta1 = dot(univec,gamma)
    delta2 = ptdf_full*gamma

    delta1_min = sum(gamma_min)
    delta1_max = sum(gamma_max)

    xmu  = 0.5*(gamma_min + gamma_max)
    xsig = 0.5*(gamma_max - gamma_min)
    delta2_min = -abs.(ptdf_full)*xsig + ptdf_full*xmu
    delta2_max = +abs.(ptdf_full)*xsig + ptdf_full*xmu

    # intermediate variables
    """ z  =    x *  y
        z1 = lambda*delta1  (scalar)
        z2 = mu1.*delta2    (vector)
        z3 = mu2.*delta2    (vector)
    """

    # auto-build the McCormick cuts
    @variable(model, z1)
    @variable(model, z2[1:length(mu1)])
    @variable(model, z3[1:length(mu2)])

    # for z1
    x  = lambda
    y  = delta1
    xl = lambda_min
    yl = delta1_min
    xu = lambda_max
    yu = delta1_max
    @constraint(model, xl*y + yl*x - xl*yl <= z1)
    @constraint(model, xu*y + yu*x - xu*yu <= z1)
    @constraint(model, z1 <= xl*y + yu*x - xl*yu)
    @constraint(model, z1 <= xu*y + yl*x - xu*yl)

    # for z2
    x  = mu1
    y  = delta2
    xl = mu1_min
    yl = delta2_min
    xu = mu1_max
    yu = delta2_max
    @constraint(model, xl.*y .+ yl.*x .- xl.*yl .<= z2)
    @constraint(model, xu.*y .+ yu.*x .- xu.*yu .<= z2)
    @constraint(model, z2 .<= xl.*y .+ yu.*x .- xl.*yu)
    @constraint(model, z2 .<= xu.*y .+ yl.*x .- xu.*yl)

    # for z3
    x  = mu2
    y  = delta2
    xl = mu2_min
    yl = delta2_min
    xu = mu2_max
    yu = delta2_max
    @constraint(model, xl.*y .+ yl.*x .- xl.*yl .<= z3)
    @constraint(model, xu.*y .+ yu.*x .- xu.*yu .<= z3)
    @constraint(model, z3 .<= xl.*y .+ yu.*x .- xl.*yu)
    @constraint(model, z3 .<= xu.*y .+ yl.*x .- xu.*yl)

    t1 = -lambda*dot(univec,pd)
    t2 = -lambda*z1

    t3 = -mu1'*ptdf_full*pd
    t4 = -sum(z2)
    t5 = -mu1'*flow_max

    t6 = mu2'*ptdf_full*pd
    t7 = sum(z3)
    t8 = -mu2'*flow_max

    t9 = -mu3'*pmax
    t10= +mu4'*pmin

    @constraint(model, 0.001 == t1+t2+t3+t4+t5+t6+t7+t8+t9+t10)
    @objective(model, Min, dot(gamma,gamma))
    optimize!(model)

    termination_status(model)
    println(value.(gamma))

    return model
end

function local_infeasibility_and_bounding(nb, ng, A, C, univec, pd, ptdf_full, flow_max, pmax, pmin, zero_perts)
    model = Model(Ipopt.Optimizer)
    @variable(model, gamma[1:nb])
    @constraint(model, gamma[1:zero_perts] .== 0.0)

    B = -univec'*(pd .+ gamma[1:nb])
    D = -[ptdf_full*(pd .+ gamma[1:nb]).+flow_max;-ptdf_full*(pd .+ gamma).+flow_max;pmax;-pmin]

    @variable(model, lambda)
    nm = size(C,1)
    @variable(model, 0.0 <= mu[1:nm])
    @constraint(model, lambda*A + mu'*C .== 0)
    @constraint(model, 0.001 ==  lambda*B + mu'*D)

    # => @constraint(model, 0.001 .<=  lambda*B + mu'*D)
    @objective(model, Min, dot(gamma,gamma))
    optimize!(model)

    println("termination status: ", termination_status(model))
    return model
end

function SDP_find_smallest_infeasibility_perturbation(nb, ng, A, univec, pd, ptdf_full, N, flow_max, pmax, pmin, zero_perts)
        
    model = Model(Mosek.Optimizer)

    nM = 1 + 1 + (2*nl+2*ng) + nb
    @variable(model, M[1:nM, 1:nM], PSD)
    @constraint(model, M[1,1]  == 1)

    # indices
    lam_idx = 2
    mu_idx  = (2+1):(2 + 2*nl+2*ng)
    gam_idx = (2 + 2*nl+2*ng + 1):(2 + 2*nl+2*ng + nb)
    mu1_idx = mu_idx[(1):(nl)]
    mu2_idx = mu_idx[(nl+1):(2*nl)]
    mu3_idx = mu_idx[(2*nl+1):(2*nl+ng)]
    mu4_idx = mu_idx[(2*nl+ng+1):(2*nl+2*ng)]

    lambda  = M[lam_idx,1]
    mu      = M[mu_idx,1]
    gamma   = M[gam_idx,1]

    mu1     = mu[(1):(nl)]
    mu2     = mu[(nl+1):(2*nl)]
    mu3     = mu[(2*nl+1):(2*nl+ng)]
    mu4     = mu[(2*nl+ng+1):(2*nl+2*ng)]

    lam_mu  = M[[lam_idx; mu_idx],[lam_idx; mu_idx]]
    gam_gam = M[gam_idx,gam_idx]
    mu_mu   = M[mu_idx,mu_idx]
    lam_gam = M[lam_idx,gam_idx]
    gam_mu1 = M[gam_idx,mu1_idx]
    gam_mu2 = M[gam_idx,mu2_idx]

    lam_mu_mu  = M[[lam_idx; mu_idx],mu_idx]
    lam_mu_gam = M[[lam_idx; mu_idx],gam_idx]
    lam_mu_lam = M[[lam_idx; mu_idx],lam_idx]

    # build the t's
    t1   = -lambda*dot(univec,pd)
    t2   = -sum(lam_gam)
    t3   = -mu1'*ptdf_full*pd
    t4   = -tr(ptdf_full*gam_mu1)
    t5   = -mu1'*flow_max
    t6   = mu2'*ptdf_full*pd
    t7   = tr(ptdf_full*gam_mu2) 
    t8   = -mu2'*flow_max
    t9   = -mu3'*pmax
    t10  = +mu4'*pmin

    # constraints
    @constraint(model, gamma[1:zero_perts] .== 0.0)
    @constraint(model, lambda*A + mu'*C    .== 0.0)
    tmat = [A' C']
    @constraint(model, tmat*lam_mu*(tmat') .== 0.0)
    @constraint(model, tmat*lam_mu_mu      .== 0.0)
    @constraint(model, tmat*lam_mu_gam     .== 0.0)
    @constraint(model, tmat*lam_mu_lam     .== 0.0)

    @constraint(model, 0.0 .<= mu)
    @constraint(model, 0.0 .<= mu_mu)

    @constraint(model, 0.001 == t1+t2+t3+t4+t5+t6+t7+t8+t9+t10 )
    @objective(model, Min, tr(gam_gam))
    optimize!(model)

    println(termination_status(model))

    return model
end

function select_slack_gen(basic_network_data)
    @warn("depricated, because we want the bus and gen index at the same time")
    # pick the generator with non-zero capacity and the smallest index
    ii = 1
    idx = 0
    for (key,val) in basic_network_data["gen"]
        if val["pmax"] > 0   
            if ii == 1
                idx = val["gen_bus"]
                ii += 1
            end
            if idx > val["gen_bus"]
                idx = val["gen_bus"]
            end
        end
    end

    return idx
end

function callback_log(cb_data, cb_where::Cint)
    MIP = 3
    if cb_where == MIP
        # println(cb_where)
        runtimeP = Ref{Cdouble}()
        objbstP  = Ref{Cdouble}()
        objbndP  = Ref{Cdouble}()
        GRBcbget(cb_data, cb_where, GRB_CB_RUNTIME, runtimeP)
        GRBcbget(cb_data, cb_where, GRB_CB_MIP_OBJBST, objbstP)
        GRBcbget(cb_data, cb_where, GRB_CB_MIP_OBJBND, objbndP)
        best  = objbstP[]
        bound = objbndP[]
        # push data
        push!(bestlog, best)
        push!(boundlog, bound)
        push!(timelog, runtimeP[])
        # don't save the gap:           
              # gap = abs((objbstP[] - objbndP[]) / objbstP[])
              # push!(gaplog, gap)
    end
    return
end

function filterlog(timelog, bestlog, boundlog)
    # filter out the "-1.0e100" values (no solution found)
    inds           = abs.(bestlog) .> 1e15
    if length(inds) > 0
        bestlog[inds] .= 0
        max_ind        = argmax(abs.(bestlog))
        bestlog[inds] .= bestlog[max_ind]
    end
    gaplog         = 100*abs.((bestlog - boundlog) ./ bestlog)

    log = Dict(
          :gaplog   => convert(Array{Float64},gaplog),
          :timelog  => convert(Array{Float64},timelog),
          :bestlog  => convert(Array{Float64},bestlog),
          :boundlog => convert(Array{Float64},boundlog))

    return log
end

function write_hdf5data(testdata_file, farkas_log)
    fid = h5open(testdata_file, "w") do file
        write(file, "gaplog",  farkas_log[:gaplog])
        write(file, "timelog", farkas_log[:timelog])
        write(file, "boundlog", farkas_log[:boundlog])
        write(file, "bestlog", farkas_log[:bestlog])
    end
end

function write_squeeze_hdf5data(testdata_file, farkas_log, t_final, delta2_final)
    fid = h5open(testdata_file, "w") do file
        write(file, "gaplog",  farkas_log[:gaplog])
        write(file, "timelog", farkas_log[:timelog])
        write(file, "boundlog", farkas_log[:boundlog])
        write(file, "bestlog", farkas_log[:bestlog])
        write(file, "t", t_final)
        write(file, "delta2", delta2_final)
    end
end

function read_hdf5data(testdata_file)
    fid      = h5open(testdata_file, "r")
    gaplog   = read(fid, "gaplog")  
    timelog  = read(fid, "timelog") 
    boundlog = read(fid, "boundlog")
    bestlog  = read(fid, "bestlog") 
    close(fid)

    return gaplog, timelog, boundlog, bestlog
end

function maximize_p0_margin(A,B,b)
    model = Model(Gurobi.Optimizer)
    npg = size(A,2)
    @variable(model, pg[1:npg])
    @variable(model, t)
    @constraint(model, A*pg + b .<= 0.0)
    for ii in 1:length(b)
        if ~all(iszero, B[ii,:])
            @constraint(model, dot(A[ii,:],pg) + b[ii] <= t)
        end
    end

    @objective(model, Min, t)
    optimize!(model)

    return value.(pg)
end

function initialize_local_control(A,B,b,pg0)
    model = Model(Gurobi.Optimizer)
    @variable(model, t)
    ab    = (A*pg0 + b).^2
    nc = size(A,1)
    vv = zeros(nc)
    for ii in 1:nc
        vv[ii] = B[ii,:]'B[ii,:]
        @constraint(model, t*vv[ii] <= ab[ii])
    end
    @objective(model, Max, t)
    optimize!(model)

    return value(t), vv
end