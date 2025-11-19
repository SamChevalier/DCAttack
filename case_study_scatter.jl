using HDF5
using Plots
using Random
using LinearAlgebra
using InvertedIndices
using PowerModels, PGLib
using JuMP, Ipopt, Gurobi

include("./src/squeeze_functions.jl")
Random.seed!(1)

# call case
case = "pglib_opf_case5_pjm.m"

# solve Farkas' lemma attack
MIPGap    = 0.001
TimeLimit = 600.00

network_data       = pglib(case)
basic_network_data = PowerModels.make_basic_network(network_data)
zero_nonlinear_costs!(basic_network_data)
basic_network_data_base = deepcopy(basic_network_data)

# locally solve Farkas' lemma
num_perts   = 1000
delta_solns = zeros(3,num_perts)
new_loads   = zeros(3,num_perts)

for ii in 1:num_perts
    if ii > 1
        basic_network_data["load"]["1"]["pd"] = copy(basic_network_data_base["load"]["1"]["pd"]*(1+0.1*randn()))
        basic_network_data["load"]["2"]["pd"] = copy(basic_network_data_base["load"]["2"]["pd"]*(1+0.1*randn()))
        basic_network_data["load"]["3"]["pd"] = copy(basic_network_data_base["load"]["3"]["pd"]*(1+0.1*randn()))

        new_loads[1,ii] = copy(basic_network_data["load"]["1"]["pd"])
        new_loads[2,ii] = copy(basic_network_data["load"]["2"]["pd"])
        new_loads[3,ii] = copy(basic_network_data["load"]["3"]["pd"])
    else
        new_loads[1,1] = copy(basic_network_data_base["load"]["1"]["pd"]);
        new_loads[2,1] = copy(basic_network_data_base["load"]["2"]["pd"]);
        new_loads[3,1] = copy(basic_network_data_base["load"]["3"]["pd"]);
    end

    pm_result, model, model_fl, sys = solve_dcopf(basic_network_data)
    A,B,c,nb = copy(sys[:A]),copy(sys[:B]),copy(sys[:c]),copy(sys[:nb])

    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "MIPGap", 0.01)
    set_optimizer_attribute(model, "TimeLimit", 1800.0)

    num_p   = size(A,2)
    num_mu  = size(A,1)
    n_perts = size(B,2)
    @variable(model, delta[1:n_perts])
    @variable(model, p[1:num_p])
    @variable(model, mu[1:num_mu])
    @constraint(model, 0.0 .<= mu .<= 1.0)
    @constraint(model, A'*mu .== 0.0)
    @constraint(model, sum(mu) == 1.0)
    @constraint(model, A*p + B*delta + c .<= 0.0001)
    @constraint(model, mu.*(A*p + B*delta + c) .== 0.0)
    @objective(model, Min, dot(delta,delta))
    optimize!(model)

    delta_solns[:,ii] .= value.(delta)
end

n  = zeros(num_perts)
nl = zeros(num_perts)
for ii in 1:num_perts
    n[ii]  = norm(delta_solns[:,ii])
    nl[ii] = norm(new_loads[:,ii])
end

# %% now, write the data
testdata_file = "./dcopf/data/CS2.h5"
fid = h5open(testdata_file, "w") do file
    write(file, "perturbation",  n)
    write(file, "load",         nl)
end

# %% now, call the data :)
fid             = h5open(testdata_file, "r")
perturbation    = read(fid, "perturbation")  
load            = read(fid, "load") 
close(fid)

using Plots.PlotMeasures
c1 = 165/256
c2 = 42/256
c3 = 42/256
redd = RGB(c1,c2,c3)
plot(load, perturbation, seriestype=:scatter, size=(600,200),xlabel="total system load (pu)", ylabel="attack size (pu)", color=redd, label="",bottom_margin=3mm,left_margin=3mm)