using HDF5
using Plots
using LinearAlgebra
using InvertedIndices
using PowerModels, PGLib
using JuMP, Ipopt, Gurobi

using Random
Random.seed!(1)

include("./src/squeeze_functions.jl")

# define test cases and make sure they are parsed correctly
cases =     ["pglib_opf_case5_pjm.m";
             "pglib_opf_case14_ieee.m";
             "pglib_opf_case24_ieee_rts.m"
             "pglib_opf_case30_as.m";
             "pglib_opf_case57_ieee.m";
             "pglib_opf_case60_c.m";
             "pglib_opf_case118_ieee.m"]

# test parsing of the test cases
test_parsing(cases)

# %% ========= defense sequence! ========= %% #
ii = 1
best_control = zeros(length(cases))
tmax_ipopt   = 600.0   
MIPGap       = 0.001
TimeLimit    = 1800.00

for case in cases

    network_data       = pglib(case)
    basic_network_data = PowerModels.make_basic_network(network_data)
    zero_nonlinear_costs!(basic_network_data)

    # solve dcopf
    pm_result, model, model_fl, sys = solve_dcopf(basic_network_data)

    # locally solve Farkas' lemma
    A,B,c,nb = copy(sys[:A]),copy(sys[:B]),copy(sys[:c]),copy(sys[:nb])

    # get locally robust generation dispatch
    pg0    = maximize_p0_margin(A, B, c)
    t0, v0 = initialize_local_control(A, B, c, pg0)
    term_stat, G0_lcl, pg0_lcl, v0_lcl, t0_lcl = solve_control_local(A, B, c, pg0, t0, v0, tmax_ipopt; random_start=false)
    if Int(term_stat) ∉ [1;4;7]
        @warn("ipopt failed!")
    end

    # initialize 
    x0 = Dict(:mu    => 0,
              :delta => 0,
              :pg0   => pg0_lcl,
              :t     => t0_lcl,
              :v     => v0_lcl,
              :G     => G0_lcl)

    # solve squeeze log
    control_log = solve_control(A, B, c, MIPGap, TimeLimit, nb, x0; init=true, extra_string="_final")

    best_control[ii] = control_log[:bestlog][end]
    ii = ii + 1
end

# %% ========= attack sequence! ========= %% #
ii = 1
best_attack  = zeros(length(cases))
tmax_ipopt   = 600.0   
MIPGap       = 0.001
TimeLimit    = 1800.00

for case in cases
    network_data       = pglib(case)
    basic_network_data = PowerModels.make_basic_network(network_data)
    zero_nonlinear_costs!(basic_network_data)

    # solve dcopf
    pm_result, model, model_fl, sys = solve_dcopf(basic_network_data)

    # locally solve Farkas' lemma
    A,B,c,nb     = copy(sys[:A]),copy(sys[:B]),copy(sys[:c]),copy(sys[:nb])
    num_mu       = size(A,1)
    n_perts      = size(B,2)

    # re-set
    Random.seed!(1)

    attack_obj = Inf
    mu_fl      = 0
    delta_fl   = 0
    for kk in 1:5
        if kk == 1
            initialize_start = false
            delta0 = 0.0
            mu0    = 0.0
        else
            initialize_start = true
            delta0 = 0.0001*randn(n_perts)
            mu0    = 10.0*rand(num_mu)
        end

        term_stat, delta_lcl, mu_lcl, obj = solve_farkas_lemma_local(A, B, c; initialize_start=initialize_start, lower_ipopt_tol=false, delta0=delta0, mu0=mu0)
        if Int(term_stat) ∉ [1;4;7]
            # decrease solver tolerance and try again
            term_stat, delta_lcl, mu_lcl, obj = solve_farkas_lemma_local(A, B, c; initialize_start=initialize_start, lower_ipopt_tol=true, delta0=delta0, mu0=mu0)
            if (Int(term_stat) ∈ [1;4;7]) && (attack_obj > obj)
                attack_obj = copy(obj)
                mu_fl      = copy(mu_lcl)
                delta_fl   = copy(delta_lcl)
            end
        else
            if attack_obj > obj
                attack_obj = copy(obj)
                mu_fl      = copy(mu_lcl)
                delta_fl   = copy(delta_lcl)
            end
        end
    end

    # initialize 
    x0 = Dict(:mu    => mu_fl,
              :delta => delta_fl,
              :pg0   => 0,
              :t     => 0,
              :v     => 0,
              :G     => 0)

    # solve farkas lemma
    faraks_log = solve_farkas_lemma(A, B, c, MIPGap, TimeLimit, nb, x0; init=true, extra_string="_final")

    best_attack[ii]  = faraks_log[:bestlog][end]
    ii = ii + 1
end