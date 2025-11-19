using Pkg
Pkg.activate(".")

# call packages
using HDF5
using Plots
using LinearAlgebra
using InvertedIndices
using PowerModels, PGLib
using JuMP, Ipopt, Gurobi

# tests normalized Farkas' lemma 

using Random
Random.seed!(1)

include("./src/squeeze_functions.jl")

# %% ========= attack sequence! ========= %% #
cases =     ["pglib_opf_case5_pjm.m";
             "pglib_opf_case14_ieee.m";
             "pglib_opf_case24_ieee_rts.m"
             "pglib_opf_case30_as.m";
             "pglib_opf_case57_ieee.m";
             "pglib_opf_case60_c.m";
             "pglib_opf_case118_ieee.m"]

case = cases[7]
network_data       = pglib(case)
basic_network_data = PowerModels.make_basic_network(network_data)
zero_nonlinear_costs!(basic_network_data)

# solve dcopf
pm_result, model, model_fl, sys = solve_dcopf(basic_network_data)

# locally solve Farkas' lemma
A,B,c,nb = copy(sys[:A]),copy(sys[:B]),copy(sys[:c]),copy(sys[:nb])

model = Model(Gurobi.Optimizer)
tg = 172800.0 - 500.0
set_optimizer_attribute(model, "MIPGap", 0.01)
set_optimizer_attribute(model, "TimeLimit", tg)
set_optimizer_attribute(model, "MIPFocus", 3)

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