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

# %% ========= attack sequence! ========= %% #

case = "pglib_opf_case118_ieee.m"
network_data       = pglib(case)
basic_network_data = PowerModels.make_basic_network(network_data)
zero_nonlinear_costs!(basic_network_data)

# solve dcopf
pm_result, model, model_fl, sys = solve_dcopf(basic_network_data)

# locally solve Farkas' lemma
A,B,b,nb = copy(sys[:A]),copy(sys[:B]),copy(sys[:b]),copy(sys[:nb])

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
@constraint(model, A*p + B*delta + b .<= epsilon)
@constraint(model, mu.*(A*p + B*delta + b) .== 0.0)

@constraint(model, -1 .<= delta .<= 1)

@objective(model, Min, dot(delta,delta))

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
extra_string = "_final_reformulation"
testdata_file = "./data/farkas_"*string(nb)*"bus"*extra_string*".h5"
write_hdf5data(testdata_file, faraks_log)

println(objective_value(model))

