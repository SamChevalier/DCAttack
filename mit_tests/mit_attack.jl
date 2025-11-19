using HDF5
using Plots
using LinearAlgebra
using InvertedIndices
using PowerModels, PGLib
using JuMP, Ipopt, Gurobi

# This script kicks off a branch-and-bound search for the
# smallest possible attack (to show that attacks are hard)

using Random
Random.seed!(1)

include("../src/squeeze_functions.jl")

# %% ========= attack sequence! ========= %% #
case = "pglib_opf_case57_ieee.m"

function run_attack(case::String)
    network_data       = pglib(case)
    basic_network_data = PowerModels.make_basic_network(network_data)

    zero_nonlinear_costs!(basic_network_data)

    # solve dcopf
    pm_result, model, model_fl, sys = solve_dcopf(basic_network_data)

    # locally solve Farkas' lemma
    A,B,c,nb = copy(sys[:A]),copy(sys[:B]),copy(sys[:c]),copy(sys[:nb])

    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "MIPGap", 0.01)
    set_optimizer_attribute(model, "TimeLimit", 100.0)

    num_mu  = size(A,1)
    n_perts = size(B,2)
    @variable(model, delta[1:n_perts])
    @variable(model, mu[1:num_mu])
    @constraint(model, 0.0 .<= mu)
    @constraint(model, A'*mu .== 0.0)
    @constraint(model, mu'*(B*delta+c) .== 0.0001)
    @constraint(model, sum(mu) == 1.0)
    @objective(model, Min, dot(delta,delta))

    # optimize
    optimize!(model)
end

run_attack(case)
