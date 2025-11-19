using HDF5
using Plots
using Random
using LinearAlgebra
using InvertedIndices
using Mosek, MosekTools
using PowerModels, PGLib
using JuMP, Ipopt, Gurobi

include("./squeeze_functions.jl")

# %% call case
case = "pglib_opf_case5_pjm.m"

# solve Farkas' lemma attack
MIPGap    = 0.001
TimeLimit = 600.00

network_data       = pglib(case)
basic_network_data = PowerModels.make_basic_network(network_data)
zero_nonlinear_costs!(basic_network_data)

# solve dcopf
pm_result, model, model_fl, sys = solve_dcopf(basic_network_data)

# %% locally solve Farkas' lemma
A,B,b,nb = copy(sys[:A]),copy(sys[:B]),copy(sys[:b]),copy(sys[:nb])
# get locally robust generation dispatch
pg0    = maximize_p0_margin(A, B, b)
t0, v0 = initialize_local_control(A, B, b, pg0)

num_mu  = size(A,1)
n_perts = size(B,2)
nx      = size(A,2)
nc      = size(A,1)
tmax_ipopt = 100.0
term_stat, G0_lcl, pg0_lcl, v0_lcl, t0_lcl = solve_control_local(A, B, b, pg0, t0, v0, tmax_ipopt; random_start=false)

println(t0_lcl)

# initialize 
Random.seed!(1)
x0 = Dict(  :mu    => randn(),
            :delta => randn(),
            :pg0   => pg0_lcl + 0*randn(nx),
            :t     => t0_lcl  + 0*randn(),
            :v     => v0_lcl  + 0*randn(nc),
            :G     => G0_lcl  + 0*randn(nx,n_perts))
extra_string = "_casestudy"

TimeLimit = 100.00
# => control_log = solve_control(A, B, b, MIPGap, TimeLimit, nb, x0; init=true, extra_string=extra_string)
# => faraks_log  = solve_farkas_lemma(A, B, b, MIPGap, TimeLimit, nb, x0; init=true, extra_string=extra_string)

# %% ===========
nb = 5

testdata_file_farkas = "./dcopf/data/farkas_"*string(nb)*"bus_casestudy.h5"
testdata_file_control = "./dcopf/data/control_"*string(nb)*"bus_casestudy.h5"
gaplog5_farkas,  timelog5_farkas,  boundlog5_farkas,  bestlog5_farkas = read_hdf5data(testdata_file_farkas)
gaplog5_control, timelog5_control, boundlog5_control, bestlog5_control = read_hdf5data(testdata_file_control)

# %% ===
gr()

scaled_upper_bound =  boundlog5_control[16:10:end]/100000
scaled_upper_bound[1:1000] = 2*scaled_upper_bound[1:1000]
scaled_upper_bound .= 20.0
plot(timelog5_control[16:10:end], bestlog5_control[16:10:end], color = :steelblue, width = 8, label = "Defense Incumbent",legendfontsize=9, linealpha = 0.75)
#plot!(timelog5_control[16:10:end], scaled_upper_bound, color = :steelblue, width = 2, linestyle = :dashdot, label = "Attack Bound (Scaled Down)")
plot!(timelog5_control[16:10:end], bestlog5_control[16:10:end],fillrange = scaled_upper_bound, fillalpha=0.2, color=:steelblue, label="Defense Gap", legend = :topright)


c1 = 165/256
c2 = 42/256
c3 = 42/256
redd = RGB(c1,c2,c3)

plot!(timelog5_farkas[3:10:end], bestlog5_farkas[3:10:end], color=redd, width = 2, label = "Attack Incumbent")
plot!(timelog5_farkas[3:10:end], boundlog5_farkas[3:10:end], color=redd, linestyle = :dashdot, width = 2, label = "Attack Bound")
plot!(timelog5_farkas[3:10:end], boundlog5_farkas[3:10:end], fillrange=bestlog5_farkas[3:10:end], fillalpha=0.2, label="Attack Gap", color=redd)

annotate!([(20, 13, ("smallest adversarial attack", 9))])
annotate!([(20, 11.8, ("proved here!", 9))])
plot!([15,1.5],[10.5,6.9],arrow=true,color=:black,linewidth=1.5,label="", xlim = [-1,100])

r = 0.5         # radius
θ = range(0, 2π, length=500)

# Parametric coordinates
x = 2*r * cos.(θ) .+ 0.2
y = r * sin.(θ) .+ 6.28
plot!(x, y, color=:black, width = 1.75, label ="",xtickfontsize = 9,ytickfontsize = 9, xlabel = "Branch and Bound solve time (seconds)", ylabel = "Bound or Incumbent")

p = plot!(size=(600,300), dpi=500)
# savefig("case_study.pdf")