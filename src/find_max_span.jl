#finding section -> maximum span length.

#using the same mix combinations.

#It's catalog all over again, but calculate back the dosage. 

#load the different L, t, Lc. 

using CSV, JSON, DataFrames

sections = JSON.parsefile(joinpath(@__DIR__,"pixel_sample.json"), dicttype = Dict{String,Float64});

all_sections =Matrix{Float64}(undef, length(sections), 3)
results = Dict{Vector, Vector}()
for i in eachindex(sections)
    val = 1000*[sections[i]["L"], sections[i]["t"], sections[i]["Lc"]]
    all_sections[i,:] = val
    results[val] = []
end

all_sections
results


push!(results[all_sections[1,:]],1)
