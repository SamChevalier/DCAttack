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

include("./src/squeeze_functions.jl")

# this script tests the bilevel reformulation

# %% ========= attack sequence! ========= %% #
cases =     ["pglib_opf_case5_pjm.m";
             "pglib_opf_case14_ieee.m";
             "pglib_opf_case24_ieee_rts.m"
             "pglib_opf_case30_as.m";
             "pglib_opf_case57_ieee.m";
             "pglib_opf_case60_c.m";
             "pglib_opf_case118_ieee.m"]
case = cases[6]
network_data       = pglib(case)
basic_network_data = PowerModels.make_basic_network(network_data)
zero_nonlinear_costs!(basic_network_data)

# solve dcopf
pm_result, model, model_fl, sys = solve_dcopf(basic_network_data)

# locally solve Farkas' lemma
A,B,c,nb = copy(sys[:A]),copy(sys[:B]),copy(sys[:c]),copy(sys[:nb])

# now, solve the reformulation ==========
model = Model(Gurobi.Optimizer)
set_optimizer_attribute(model, "MIPGap", 0.01)
tg = 172800.0 - 500.0
set_optimizer_attribute(model, "TimeLimit", tg)

num_p   = size(A,2)
num_mu  = size(A,1)
n_perts = size(B,2)
@variable(model, delta[1:n_perts])
@variable(model, p[1:num_p])

epsilon = 0.001

@variable(model, mu[1:num_mu])
@constraint(model, 0.0 .<= mu .<= 1.0)
@constraint(model, A'*mu .== 0.0)
@constraint(model, sum(mu) == 1.0)
@constraint(model, A*p + B*delta + c .<= epsilon)
@constraint(model, mu.*(A*p + B*delta + c) .== 0.0)
@objective(model, Min, dot(delta,delta))

# optimize
optimize!(model)
println(objective_value(model))