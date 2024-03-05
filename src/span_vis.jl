using JSON, DataFrames
using Makie, GLMakie

catalog = JSON.parsefile(joinpath(@__DIR__,"Results/04_03_latin_span.json"), dicttype = Dict{String,Any})

f1 = Figure(size = (1000,1000)) 
ax1 = Axis(f1[1,1], title = "Span Moment vs Shear", xlabel = "Moment controlled span [m]", ylabel = "Shear controlled span [m]") 

moment = Vector{Float64}()
shear = Vector{Float64}() 
set_color =Vector{Symbol}()
for k in keys(catalog)
    sub_catalog = catalog[k]
    for i in eachindex(sub_catalog)
        push!(moment, sub_catalog[i]["span_mu"])
        push!(shear, sub_catalog[i]["span_vu"])
        color = sub_catalog[i]["control"] == "Moment control" ? :green : :red 
        push!(set_color, color)
    end
end

scatter!(ax1, moment, shear, color = set_color)
scatter!(ax1, 0:1:30, 0:1:30)

f1


output = Vector{Dict{String, Any}}(undef, length(catalog))
for k in keys(catalog)
    sub_catalog = catalog[k]
    max_l = -1
    idx = -1
    for i in eachindex(sub_catalog)
        sub_l = sub_catalog[i]["span"]
        control_case = sub_catalog[i]["control"]
        max_l = max_l > sub_l ? max_l : sub_l
        idx = max_l > sub_l ? idx : i
    end
    output[parse(Int64, k)] = sub_catalog[idx]
end

open("src//Results//04_03_span_max_latin.json","w") do f
    JSON.print(f, output)
end

CSV.write("src//Results//04_03_span_max_latin.csv", DataFrame(output))
