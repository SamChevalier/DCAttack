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

# %% test
include("./src/squeeze_functions.jl")

MIPGap       = 0.001
TimeLimit    = 15.00

case = cases[4]
network_data       = pglib(case)
basic_network_data = PowerModels.make_basic_network(network_data)
zero_nonlinear_costs!(basic_network_data)

# solve dcopf
pm_result, model, model_fl, sys = solve_dcopf(basic_network_data)

# locally solve Farkas' lemma
A,B,b,nb     = copy(sys[:A]),copy(sys[:B]),copy(sys[:b]),copy(sys[:nb])
num_mu       = size(A,1)
n_perts      = size(B,2)

# initialize 
x0   = Dict(:mu    => 0,
            :delta => 0,
            :pg0   => 0,
            :t     => 0,
            :v     => 0,
            :G     => 0)

# solve farkas lemma
solve_farkas_lemma_record_log(A, B, b, MIPGap, TimeLimit, nb, x0; init=false)


