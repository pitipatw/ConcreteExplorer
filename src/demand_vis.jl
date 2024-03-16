using Makie, GLMakie
using DataFrames, JSON

loaded_demands = JSON.parsefile(joinpath("src//test_demands.json"), dicttype =Dict{String, Vector{Union{Int64,Float64, String}}} );

figure1 = Figure(size = (1000,1000))
ax1 = Axis(figure1[1,1])
ax2 = Axis(figure1[2,1])
scatter!(ax1, Int.(loaded_demands["idx"][1:100]), Int.(loaded_demands["mu"][1:100]))
scatter!(ax2, Int.(loaded_demands["idx"][1:100]), Int.(loaded_demands["vu"][1:100]))
figure1

loaded_demands_df = DataFrame(loaded_demands)

