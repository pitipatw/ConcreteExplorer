using Makie, GLMakie
using DataFrames, JSON

loaded_demands = JSON.parsefile(joinpath("src//test_demands.json"), dicttype =Dict{String, Vector{Union{Int64,Float64, String}}} );

figure1 = Figure(size = (1000,1000))
ax1 = Axis(figure1[1,1])
ax2 = Axis(figure1[2,1])
scatter!(ax1, Int.(loaded_demands["idx"][13:24]), Int.(loaded_demands["mu"][13:24]))
scatter!(ax2, Int.(loaded_demands["idx"][13:24]), Int.(loaded_demands["vu"][13:24]))