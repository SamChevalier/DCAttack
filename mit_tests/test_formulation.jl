using Pkg
Pkg.activate(".")

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
case               = "pglib_opf_case57_ieee.m"
network_data       = pglib(case)
basic_network_data = PowerModels.make_basic_network(network_data)
zero_nonlinear_costs!(basic_network_data)

# solve dcopf
pm_result, model, model_fl, sys = solve_dcopf(basic_network_data)

# parameterize the problem
A,B,c,nb = copy(sys[:A]),copy(sys[:B]),copy(sys[:c]),copy(sys[:nb])

testdata_file = "./mit_tests/ABc_data.h5"
fid = h5open(testdata_file, "w") do file
    write(file, "A",  A)
    write(file, "B",  B)
    write(file, "c",  c)
end

function find_smallest_attack(A, B, c; solve_strategy=:local, tmax=50.0)
    if solve_strategy == :local
        model = Model(Ipopt.Optimizer)
    elseif solve_strategy == :global
        model = Model(Gurobi.Optimizer)
        set_optimizer_attribute(model, "TimeLimit", tmax)
    else
        @warn("solve strategy not recognized!!")
    end

    @variable(model, delta[1:size(B,2)])
    @variable(model, mu[1:size(A,1)])
    @constraint(model, 0.0 .<= mu)
    @constraint(model, A'*mu .== 0.0)
    @constraint(model, mu'*(B*delta+c) .== 0.0001)
    @constraint(model, sum(mu) == 1.0)
    @objective(model, Min, dot(delta,delta))
    optimize!(model)

    return objective_value(model)
end

# %% solve for the smallest attacks
attack_size_local = find_smallest_attack(A, B, c; solve_strategy = :local)
attack_size_global = find_smallest_attack(A, B, c; solve_strategy = :global)

@info("best local attack size: $attack_size_local")
@info("best global attack size (found in given time limit): $attack_size_global")
@info("true solution (see paper): 0.0547")

# %% initialization of pg0, from our paper.
# See eq. (29) in our paper: https://arxiv.org/pdf/2507.07850

# This is similar to the Chebyshev center proposed by Dirk.

# We find the nominal generation dispatch which is (1) feasible and
# maximizes the inequality constraint margins (l-inf norm).

# We do this to warm-start in the nonlinear program we solve 
# to find the control policy which maximizes distance to infeasibility.

# note: there are only 3 non-slack generators.
pg0 = maximize_p0_margin(A, B, c)
