using Pkg
Pkg.activate(".")

# call packages
using HDF5
using Plots
using LinearAlgebra
using InvertedIndices
using PowerModels, PGLib
using JuMP, Ipopt, Gurobi

using Random
Random.seed!(1)

# %% test
include("./src/squeeze_functions.jl")

function run_118_test(TimeLimit::Float64)
    MIPGap    = 0.001
    # => TimeLimit = 100.0

    # %% ========= attack sequence! ========= %% #
    tmax_ipopt  = 600.0   

    case = "pglib_opf_case118_ieee.m"
    network_data       = pglib(case)
    basic_network_data = PowerModels.make_basic_network(network_data)
    zero_nonlinear_costs!(basic_network_data)

    # solve dcopf
    pm_result, model, model_fl, sys = solve_dcopf(basic_network_data)

    # locally solve Farkas' lemma
    A,B,b,nb     = copy(sys[:A]),copy(sys[:B]),copy(sys[:b]),copy(sys[:nb])
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

        term_stat, delta_lcl, mu_lcl, obj = solve_farkas_lemma_local(A, B, b; initialize_start=initialize_start, lower_ipopt_tol=false, delta0=delta0, mu0=mu0)
        if Int(term_stat) ∉ [1;4;7]
            # decrease solver tolerance and try again
            term_stat, delta_lcl, mu_lcl, obj = solve_farkas_lemma_local(A, B, b; initialize_start=initialize_start, lower_ipopt_tol=true, delta0=delta0, mu0=mu0)
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
    faraks_log = solve_farkas_lemma(A, B, b, MIPGap, TimeLimit, nb, x0; init=true, extra_string="_final_7day")
end

# run
tg = 172800.0 - 500.0
run_118_test(tg)

# %% ===
#using HDF5
#nb = 118
#extra_string="_final_24h"
#testdata_file = "./data/farkas_"*string(nb)*"bus"*extra_string*".h5"
#fid      = h5open(testdata_file, "r")
#gaplog   = read(fid, "gaplog")  
#timelog  = read(fid, "timelog") 
#boundlog = read(fid, "boundlog")
#bestlog  = read(fid, "bestlog") 
#close(fid)

# squeue -u $USER
# tail -f solve_118_1791051.out
# tail -n +1 -f solve_118_1791051.out
# tail -n +1 -f solve_118_1792140.out

#fid      = h5open(testdata_file, "r")
#gaplog   = read(fid, "gaplog")  
#timelog  = read(fid, "timelog") 
#boundlog = read(fid, "boundlog")
#bestlog  = read(fid, "bestlog") 
#close(fid)

#gr(show=false)
#plot(timelog/3600, bestlog)
#plot!(timelog/3600,boundlog)
#savefig("comp_h.png")
