#finding section -> maximum span length.

#using the same mix combinations.

#It's catalog all over again, but calculate back the dosage. 

#load the different L, t, Lc. 

using CSV, JSON, DataFrames, JSON
using Makie, GLMakie
include("Functions//Geometry//pixelgeo.jl")
include("Functions//capacities.jl")
include("Functions/generalfunctions.jl")
sections = JSON.parsefile(joinpath(@__DIR__,"pixel_sample.json"), dicttype = Dict{String,Float64});

all_sections =Matrix{Float64}(undef, length(sections), 3)
# results = Dict{Vector, Vector}()
for i in eachindex(sections)
    val = 1000*[sections[i]["L"], sections[i]["t"], sections[i]["Lc"]]
    all_sections[i,:] = val
    # results[val] = []
end

all_sections
ns = size(all_sections)[1]

f_all_sections = Figure(size = (500,500))
ax1 = Axis(f_all_sections[1,1], title = "L vs t", xlabel = "L [mm]", ylabel = "t [mm]") 
ax2 = Axis(f_all_sections[1,2], title = "Lc vs t",  xlabel = "Lc [mm]", ylabel = "t [mm]") 
ax3 = Axis(f_all_sections[2,1], title = "L` vs Lc", xlabel = "L [mm]", ylabel = "Lc [mm]") 
ax4 = f_all_sections[2,2] = GridLayout()

s1 = scatter!(ax1, all_sections[:,1], all_sections[:,2])
s2 = scatter!(ax2, all_sections[:,3], all_sections[:,2])
s3 = scatter!(ax3, all_sections[:,1], all_sections[:,3])
f_all_sections
#=============================================================================#
#=============================================================================#
#=============================================================================#

#load definition [kN/m]
#Concrete density = 2400 kg/m3 -> 24000 N/m3 -> 24 kN/m3
#dead load = 24*section area + 24*bay*0.1 (10 cm thick floor) kN/m
#live load = 2.4 kN/m2 - > 2.4 * bay kN/m

fc′_path = "src//Tables//fc_fiber.csv"
fc_fiber = CSV.read(fc′_path, DataFrame)
#These come in a set of (fc′,fR1, fR3)
set_fc′= convert(Array{Float64}, fc_fiber[!, :strength]) #some numbers can be Int.
set_fR1 = fc_fiber[!, :fR1]
set_fR3 = fc_fiber[!, :fR3]
set_dosage = fc_fiber[!, :dosage]

# set_fc′ = [40.0, 60.0, 80.0]
# set_fR1 = [2.97, 3.84, 5.11]
# set_fR3 = [3.50, 4.45, 5.73]

set_as = [99.0, 140.0]
set_fpe = (0:0.1:1)*0.7*1860
set_spans = 4000:500:15000

n = prod([length.([set_fc′, set_as, set_fpe, set_spans]);ns])
output = Dict{Int64,Vector{Dict{String, Union{Real,String}}}}()

global n_output = 0
for i_section in 1:ns
    L,t,Lc = all_sections[i_section,:]
    section = make_Y_layup_section(L,t,Lc)

    #viz section (ongoing)
    # for s in section.solids
    #     a = Axis(ax4[x,y])
    #     hidespines!(a)
    #     poly!(s.points)
    # end

    set_dps = L/3:L/6:L


    for i_fc′ in range(1,3)
        fc′ = set_fc′[i_fc′]
        fR1 = set_fR1[i_fc′]
        fR3 = set_fR3[i_fc′]

        for i_dps in eachindex(set_dps)
            dps = set_dps[i_dps]
            pt_pos = [-L/2 -dps ; L/2 -dps]
            for i_as in eachindex(set_as)
                pt_area = set_as[i_as]*[1.0,1.0]
                for i_fpe in eachindex(set_fpe)
                    fpe = set_fpe[i_fpe] 
                    pt_force = fpe*pt_area

                    pixelframesection = PixelFrameSection(fc′, fR1, fR3, pt_area, pt_force, pt_pos)
                    
                    #mu = get_Mu(pixelframesection)
                    #vu = get_Vu(pixelframesection)
                    as = sum(pt_area)
                    mu = get_Mu(section, fc′, as, fpe, dps, L) #kNm
                    #turn mu to kNmm
                    mu = mu*1000
                    vu = get_Vu(section, fc′, fR1, fR3, as, fpe, L,) #kN

                    #from the available mu and vu, back calculate the maximum span 

                    DL = section.area * 24 + 4000 * 0.1* .24
                    LL = 2.4*4000 #4000 is the width of the bay.

                    #1.2*DL  + 1.6 * LL  = load per m -> w 
                    w = 1.2* DL  +  1.6*LL 
                    #maximum moment in terms of w 
                    # Mu = wl^2/8
                    #Vu = wl/2 

                    l_mu = sqrt(8*mu/w)
                    l_vu = 2*vu/w 

                    span = minimum(l_mu, l)_vu 
                    status = l_mu < l_vu ? "Moment control" : "Shear control"
                        
                    this_result = Dict("fc′"=> fc′, "fR1"=> fR1, "fR3"=> fR3, 
                                    "fpe" => fpe, "as" => as, "dps" => dps,
                                    "span" => span,
                                    "L" => L , "t" => t , "Lc" => Lc,
                                    "control" => status
                                    )
                    if haskey(output, i_section)
                        push!(output[i_section],this_result)
                    else 
                        output[i_section] = [this_result]
                    end

                    end
                end
            end
        end
    end
end
println("There are ",n_output)

open("src//Results//27_02_span.json","w") do f
    JSON.print(f, output)
end

# results

# Output -> [L,t,Lc] 
# output_matrix = 
# push!(results[all_sections[1,:]],1)
