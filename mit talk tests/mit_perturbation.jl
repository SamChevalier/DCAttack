using HDF5
using Plots
using Random
using LinearAlgebra
using InvertedIndices
using Mosek, MosekTools
using PowerModels, PGLib
using JuMP, Ipopt, Gurobi

# This script tests local attacks: given a load
# perturbation, PowerModels tests for feasibility

include("../src/squeeze_functions.jl")
Random.seed!(1)

# call case
case = "pglib_opf_case57_ieee.m"

# solve Farkas' lemma attack
MIPGap    = 0.001
TimeLimit = 600.00

network_data       = pglib(case)
basic_network_data = PowerModels.make_basic_network(network_data)
zero_nonlinear_costs!(basic_network_data)
basic_network_data_base = deepcopy(basic_network_data)


# %% ======
p = plot(xlim = [0,100], ylim = [0,10], ylabel = "load perturbation 2-norm (pu)", xlabel = "trial number")
trial = 1
smallest_inf = 100.0
for trial in 1:100
    basic_network_data["load"] = deepcopy(basic_network_data_base["load"])

    kk = 1
    var = rand()
    delta = var*randn(nb)
    for (load_idx,load) in basic_network_data["load"]
        load["pd"] += copy(delta[kk])
        kk+=1
    end

    pm_result = PowerModels.solve_dc_opf(basic_network_data, Gurobi.Optimizer)
    println(pm_result["termination_status"])

    if Int(pm_result["termination_status"]) == 1
        annotate!(trial, norm(delta), text("✓", :green, 20))
    else
        annotate!(trial, norm(delta), text("✗", :red, 20))
        if smallest_inf > norm(delta)
            smallest_inf = norm(delta)
        end
    end

    display(p)
end

plot!(collect(1:100), smallest_inf*ones(100), color=:red,linestyle=:dash,)
