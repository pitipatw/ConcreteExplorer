# include("beamcalc.jl")

using CSV, DataFrames
using Makie, GLMakie
#loop through dataset.

dataset_Broyles = CSV.read("src//Tables//Compiled Concrete EPD Data Revised 11_06_2023.csv", DataFrame)
#remove one with too high Embodied Carbon (mroe than 2000)
ec_Broyles = dataset_Broyles[!, "A1-A3 Global Warming Potential (kg CO2-eq)"]

filter = ec_Broyles .> 800
deleteat!(dataset_Broyles, filter)

fc′_Broyles_psi = dataset_Broyles[!,"Concrete Compressive Strength (psi)"]
ec_Broyles = dataset_Broyles[!, "A1-A3 Global Warming Potential (kg CO2-eq)"]
#convert to MPa
fc′_Broyles_MPa = fc′_Broyles_psi*0.00689476





fig_Broyles = Figure(size = (400,400))
ax_Broyles = Axis(fig_Broyles[1,1]) 
# boxplot!(ax_hist, label, embodied_carbon_all)

#loop indices
embodied_carbon_all = zeros(length(fc′_Broyles_MPa ))
valid = zeros(length(fc′_Broyles_MPa ))
for i in eachindex(fc′_Broyles_MPa)
    fc′ = fc′_Broyles_MPa[i]
    ec_concrete = ec_Broyles[i]
    rc_section, serviceability = beam_design(fc′, ec_concrete*10^(-9));
    if isnothing(rc_section)
        embodied_carbon_all[i] = 0.0
    else
        embodied_carbon_all[i] = rc_section.embodied_carbon
    end
    valid[i] = serviceability
end

ax_beams = Axis(fig_Broyles[2,1]) 
scatter!(ax_beams, fc′_Broyles_MPa,embodied_carbon_all)


#define groups as done in the paper
label = Vector{Int64}(undef, length(fc′_Broyles_MPa))
for i in eachindex(fc′_Broyles_MPa)
    fc′ = fc′_Broyles_psi[i]
    if fc′ <2000
        label[i] = 1
    elseif fc′< 3000
        label[i] = 2
    elseif fc′< 4000
        label[i] = 3
    elseif fc′< 5000
        label[i] = 4
    elseif fc′< 6000
        label[i] = 5
    elseif fc′< 8000
        label[i] = 6
    elseif fc′< 10000
        label[i] = 7
    elseif fc′< 12000
        label[i] = 8
    else
        label[i] = 9
    end
end


boxplot!(ax_Broyles,label,ec_Broyles, color = :red)
ax_hist = Axis(fig_Broyles[1,2])
ax_hist.xticks=  1:9
ax_hist.xtickformat =  x -> ["100","2","3","4","5","6","7","8","9"][Int.(x)]
boxplot!(ax_hist, label, embodied_carbon_all)
# boxplot!(ax_hist, [3,1], [1,2], label = ["cat", "dog"])

fig_Broyles