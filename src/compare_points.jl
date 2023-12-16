include("beamcalc.jl")

using CSV, DataFrames
#loop through dataset.

dataset_jon = CSV.read("src/Compiled Concrete EPD Data Revised 11_06_2023.csv", DataFrame)
fc′_jon_psi = dataset_jon[!,"Concrete Compressive Strength (psi)"]
#convert to MPa
fc′_jon_MPa = fc′_jon_psi*0.00689476
ec_jon = dataset_jon[!, "A1-A3 Global Warming Potential (kg CO2-eq)"]

f_jon = Figure(size = (400,400))
ax_jon = Axis(f_jon[1,1]) 
scatter!(ax_jon, fc′_jon_MPa,ec_jon)


#loop indices
embodied_carbon_all = zeros(length(fc′_jon_MPa ))
valid = zeros(length(fc′_jon_MPa ))
for i in eachindex(fc′_jon_MPa)
    fc′ = fc′_jon_MPa[i]
    ec_concrete = ec_jon[i]
    rc_section, serviceability = beam_design(fc′, ec_concrete*10^(-9));
    if isnothing(rc_section)
        embodied_carbon_all[i] = 0.0
    else
        embodied_carbon_all[i] = rc_section.embodied_carbon
    end
    valid[i] = serviceability
end

# f_beams = Figure(size = (400,400))
ax_beams = Axis(f_jon[1,2]) 
scatter!(ax_beams, fc′_jon_MPa,embodied_carbon_all)


#define groups as done in the paper
label = Vector{Int64}(undef, length(fc′_jon))
for i in eachindex(fc′_jon_MPa)
    fc′ = fc′_jon_psi[i]
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

ax_hist = Axis(f_jon[2,2])
boxplot(label, embodied_carbon_all)
