using HDF5
using Plots
using Random
using LinearAlgebra
using InvertedIndices
using PowerModels, PGLib
using JuMP, Ipopt, Gurobi

# include!
include("./src/squeeze_functions.jl")

nb = 5
testdata_file_farkas = "./data/farkas_"*string(nb)*"bus_final.h5"
testdata_file_control = "./data/control_"*string(nb)*"bus_final.h5"
gaplog5_farkas,  timelog5_farkas,  boundlog5_farkas,  bestlog5_farkas = read_hdf5data(testdata_file_farkas)
gaplog5_control, timelog5_control, boundlog5_control, bestlog5_control = read_hdf5data(testdata_file_control)

nb = 14
testdata_file_farkas = "./data/farkas_"*string(nb)*"bus_final.h5"
testdata_file_control = "./data/control_"*string(nb)*"bus_final.h5"
gaplog14_farkas,  timelog14_farkas,  boundlog14_farkas,  bestlog14_farkas = read_hdf5data(testdata_file_farkas)
gaplog14_control, timelog14_control, boundlog14_control, bestlog14_control = read_hdf5data(testdata_file_control)

nb = 24
testdata_file_farkas = "./data/farkas_"*string(nb)*"bus_final.h5"
testdata_file_control = "./data/control_"*string(nb)*"bus_final.h5"
gaplog24_farkas,  timelog24_farkas,  boundlog24_farkas,  bestlog24_farkas = read_hdf5data(testdata_file_farkas)
gaplog24_control, timelog24_control, boundlog24_control, bestlog24_control = read_hdf5data(testdata_file_control)

nb = 30
testdata_file_farkas = "./data/farkas_"*string(nb)*"bus_final.h5"
testdata_file_control = "./data/control_"*string(nb)*"bus_final.h5"
gaplog30_farkas,  timelog30_farkas,  boundlog30_farkas,  bestlog30_farkas = read_hdf5data(testdata_file_farkas)
gaplog30_control, timelog30_control, boundlog30_control, bestlog30_control = read_hdf5data(testdata_file_control)

nb = 57
testdata_file_farkas = "./data/farkas_"*string(nb)*"bus_final.h5"
testdata_file_control = "./data/control_"*string(nb)*"bus_final.h5"
gaplog57_farkas,  timelog57_farkas,  boundlog57_farkas,  bestlog57_farkas = read_hdf5data(testdata_file_farkas)
gaplog57_control, timelog57_control, boundlog57_control, bestlog57_control = read_hdf5data(testdata_file_control)

nb = 60
testdata_file_farkas = "./data/farkas_"*string(nb)*"bus_final.h5"
testdata_file_control = "./data/control_"*string(nb)*"bus_final.h5"
gaplog60_farkas,  timelog60_farkas,  boundlog60_farkas,  bestlog60_farkas = read_hdf5data(testdata_file_farkas)
gaplog60_control, timelog60_control, boundlog60_control, bestlog60_control = read_hdf5data(testdata_file_control)

nb = 118
testdata_file_farkas = "./data/farkas_"*string(nb)*"bus_final.h5"
testdata_file_control = "./data/control_"*string(nb)*"bus_final.h5"
gaplog118_farkas,  timelog118_farkas,  boundlog118_farkas,  bestlog118_farkas = read_hdf5data(testdata_file_farkas)
gaplog118_control, timelog118_control, boundlog118_control, bestlog118_control = read_hdf5data(testdata_file_control)

# %% build the table
println("column 1")
println(boundlog5_farkas[end])
println(boundlog14_farkas[end])
println(boundlog24_farkas[end])
println(boundlog30_farkas[end])
println(boundlog57_farkas[end])
println(boundlog60_farkas[end])
println(boundlog118_farkas[end])
println()

println("column 2")
println(bestlog5_farkas[end])
println(bestlog14_farkas[end])
println(bestlog24_farkas[end])
println(bestlog30_farkas[end])
println(bestlog57_farkas[end])
println(bestlog60_farkas[end])
println(bestlog118_farkas[end])
println()

println("column 3")
println(bestlog5_control[end])
println(bestlog14_control[end])
println(bestlog24_control[end])
println(bestlog30_control[end])
println(bestlog57_control[end])
println(bestlog60_control[end])
println(bestlog118_control[end])
println()

println("column 4")
println(boundlog5_control[end])
println(boundlog14_control[end])
println(boundlog24_control[end])
println(boundlog30_control[end])
println(boundlog57_control[end])
println(boundlog60_control[end])
println(boundlog118_control[end])
println()

println("column 5")
println(bestlog5_farkas[1])
println(timelog5_farkas[1])
println(bestlog5_control[1])
println(timelog5_control[1])
println()

println(bestlog14_farkas[1])
println(timelog14_farkas[1])
println(bestlog14_control[1])
println(timelog14_control[1])
println()

println(bestlog24_farkas[1])
println(timelog24_farkas[1])
println(bestlog24_control[1])
println(timelog24_control[1])
println()

println(bestlog30_farkas[1])
println(timelog30_farkas[1])
println(bestlog30_control[1])
println(timelog30_control[1])
println()

println(bestlog57_farkas[1])
println(timelog57_farkas[1])
println(bestlog57_control[1])
println(timelog57_control[1])
println()

println(bestlog60_farkas[7])
println(timelog60_farkas[7])
println(bestlog60_control[1])
println(timelog60_control[1])
println()

println(bestlog118_farkas[1])
println(timelog118_farkas[1])
println(bestlog118_control[1])
println(timelog118_control[1])

# %% =============
using PGLib
using PowerModels
using PowerPlots

case = "pglib_opf_case57_ieee.m"
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
Random.seed!(3)
delta0 = randn(n_perts)
mu0    = rand(num_mu)

term_stat, delta_lcl, mu_lcl, obj = solve_farkas_lemma_local(A, B, b; initialize_start=true, lower_ipopt_tol=false, delta0=delta0, mu0=mu0)
println(obj)

percent_load_change = 100*(sys[:N_delta]*delta_lcl)/sum(sys[:pd])

# %% get plotting dimensions
data       = deepcopy(network_data)
data       = layout_network(data; layout_algorithm=kamada_kawai, fixed=false)
max_change = maximum(percent_load_change)
min_change = minimum(percent_load_change)
normalize_c(x) = ((x - meanv)/diffv + 1)/2

c1 = 165/256
c2 = 42/256
c3 = 42/256
redd = RGB(c1,c2,c3)
meanv = 0.5*(max_change + min_change)
diffv = 0.5*(max_change - min_change)

cmap = cgrad([redd, :white, :green], [0, normalize_c(0), 1.0])
bus_idx = []
x = []
y = []
cc = []
x_sorted = zeros(sys[:nb])
y_sorted = zeros(sys[:nb])

for (key,val) in data["bus"]
    idx = val["bus_i"]
    push!(bus_idx, idx)
    push!(y,val["ycoord_1"])
    push!(x,val["xcoord_1"])
    push!(cc,percent_load_change[idx])
    x_sorted[idx] = val["xcoord_1"]
    y_sorted[idx] = val["ycoord_1"]
end

# %% ===
p = plot()
for (key,val) in data["branch"]
    to_bus = val["t_bus"]
    fr_bus = val["f_bus"]
    x1 = x_sorted[to_bus]
    y1 = y_sorted[to_bus]
    x2 = x_sorted[fr_bus]
    y2 = y_sorted[fr_bus]
    plot!([x1,x2], [y1,y2], label = "", color=:black, width = 2.5, linealpha = 0.5)
    display(p)
end

scatter!(x,y; zcolor=cc, color=cmap, markerstrokecolor = :black,markerstrokewidth = 2, axis=false, framestyle=:none, colorbar_title = "Load change (% of total load)", ms = 8, label = "", colorbar_titlefontsize = 10)
# savefig("heat3.pdf")
