using Makie, GLMakie, CairoMakie
using JSON
using DataFrames, CSV

"""
Visualizing the design result
"""

demands = 
designs = JSON.parsefile(joinpath(@__DIR__,"Results/designs_results_13_01.json"), dicttype = Dict{String,Vector{Vector{Float64}}});

ne = length(designs)
println("There are $ne elements.")

fig = Figure(size = (1200,1200))

axs = [Axis(fig[i,j] for i in 1:10 for j in 1:10)]

