#finding section -> maximum span length.

#using the same mix combinations.

#It's catalog all over again, but calculate back the dosage. 

#load the different L, t, Lc. 

using CSV, JSON, DataFrames
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
ax1 = Axis(f_all_sections[1,1], title = "L vs Lc", xlabel = "L [mm]", ylabel = "Lc [mm]") 
ax2 = Axis(f_all_sections[1,2], title = "t vs Lc",  xlabel = "t [mm]", ylabel = "L [mm]") 
ax3 = Axis(f_all_sections[2,1], title = "L` vs t", xlabel = "L [mm]", ylabel = "t [mm]") 
ax4 = f_all_sections[2,2] = GridLayout()

s1 = scatter!(ax1, all_sections[:,1], all_sections[:,2])
s2 = scatter!(ax2, all_sections[:,3], all_sections[:,2])
s3 = scatter!(ax3, all_sections[:,1], all_sections[:,3])
#=============================================================================#
#=============================================================================#
#=============================================================================#

set_fc′ = [40.0, 60.0, 80.0]
set_fR1 = [2.97, 3.84, 5.11]
set_fR3 = [3.50, 4.45, 5.73]

set_as = 2*[99.0, 140.0]
set_fpe = 0:1860/2:1860
set_spans = 1000:5000:12000


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
                    vu = get_Vu(section, fc′, fR1, fR3, as, fpe, L,) #kN

                    #turn mu to kNmm
                    mu = mu*1000
                    #Get max span length from point load.
                    for i_spasn in eachindex(set_spans) 
                        global n_output +=1 

                        span = set_spans[i_span]
                        #Moment = PL/4
                        max_f_mu = 4*mu/span

                        #Shear = P/2
                        max_f_vu = 2*vu

                        check = max_f_mu > max_f_vu 
                        control = check ? "Shear control" : "Moment control"
                        load = check ? max_f_vu : max_f_mu
                        
                        this_result = Dict("fc′"=> fc′, "fR1"=> fR1, "fR3"=> fR3, 
                                        "fpe" => fpe, "as" => as, "dps" => dps,
                                        "load" => load, "span" => span,
                                        "L" => L , "t" => "t" , "Lc" => Lc,
                                        "control" => control
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


#visualziaing the result. 

f_output = Figure(size = (500,500))
ax1 = Axis(f_output[1,1], title = "ID vs Sections")
ax2 = Axis(f_output[1,2], title = "Span vs Load")
loads = Vector{Float64}(undef, n_output)
idx = Vector{Float64}(undef, n_output)
span = Vector{Symbol}(undef, n_output)
check = Vector{Symbol}(undef, n_output)
global count = 0
for i in 1:100
    # L,t,Lc = all_sections[i,:]
    for configs in output[i]
        global count +=1 
        loads[count] = configs["load"]
        stat  = configs["control"]
        if stat == "Shear control"
            check[count] = :red
        elseif stat == "Moment control"
            check[count] = :blue 
        end

        idx[count] = i


    end
end

scatter!(ax1, idx, loads, color = check)
scatter!()
f_output


f_all_sections

# find capacity of the section 








# results

# Output -> [L,t,Lc] 
# output_matrix = 
# push!(results[all_sections[1,:]],1)
