#finding section -> maximum span length.

#using the same mix combinations.

#It's catalog all over again, but calculate back the dosage. 

#load the different L, t, Lc. 

using CSV, JSON, DataFrames

sections = JSON.parsefile(joinpath(@__DIR__,"pixel_sample.json"), dicttype = Dict{String,Float64});

all_sections =Matrix{Float64}(undef, length(sections), 3)

for i in eachindex(sections)
    all_sections[i,:] = 1000*[sections[i]["L"], sections[i]["t"], sections[i]["Lc"]]
end

all_sections