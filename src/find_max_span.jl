#finding section -> maximum span length.

#using the same mix combinations.

#It's catalog all over again, but calculate back the dosage. 

#load the different L, t, Lc. 

using CSV, JSON, DataFrames
using Makie, GLMakie
include("Functions//Geometry//pixelgeo.jl")
include("Functions//capacities.jl")
sections = JSON.parsefile(joinpath(@__DIR__,"pixel_sample.json"), dicttype = Dict{String,Float64});

all_sections =Matrix{Float64}(undef, length(sections), 3)
results = Dict{Vector, Vector}()
for i in eachindex(sections)
    val = 1000*[sections[i]["L"], sections[i]["t"], sections[i]["Lc"]]
    all_sections[i,:] = val
    results[val] = []
end

all_sections
ns = size(all_sections)[1]

f_all_sections = Figure(size = (500,500))
ax1 = Axis(f_all_sections[1,1], title = "L vs Lc") 
ax2 = Axis(f_all_sections[1,2], title = "L vs t") 
ax3 = Axis(f_all_sections[2,1], title = "t vs Lc") 
ax4 = f_all_sections[2,2] = GridLayout()

s1 = scatter!(ax1, all_sections[:,1], all_sections[:,2])
s2 = scatter!(ax2, all_sections[:,1], all_sections[:,3])
s3 = scatter!(ax3, all_sections[:,2], all_sections[:,3])
xgrid = ns/10 
ygrid = ns/10

set_fc′ = [40.0, 60.0, 80.0]
set_fR1 = [2.97, 3.84, 5.11]
set_fR3 = [3.50, 4.45, 5.73]

set_as = [99.0, 140.0]
set_fpe = 0:1860/10:1860

global count = 0 
for x in range(1,10)
    for y in range(1,10)
        global count += 1 
        span = 0
        L,t,Lc = all_sections[count,:]
        section = make_Y_layup_section(L,t,Lc)

        #viz section (ongoing)
        # for s in section.solids
        #     a = Axis(ax4[x,y])
        #     hidespines!(a)
        #     poly!(s.points)
        # end

        set_dps = 0:L/10:L
  

        for i_fc′ in range(1,3):
            fc′[i] = set_fc′
            fR1[i] = set_fR1
            fR3[i] = set_fR3

            for i_dps in eachindex(set_dps)
                dps = set_dps[i]
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
                        mu = get_Mu(section, fc′, as, fpe, dps, L)
                        vu = get_Vu(section, fc′, fR1, fR3, as, fpe, L,)

                        #Get max span length from point load.

                        
                    end
                end
            end
        end
    end



            #create a PixelFrame section

            
            

    end
end



f_all_sections

# find capacity of the section 








# results

# Output -> [L,t,Lc] 
# output_matrix = 
# push!(results[all_sections[1,:]],1)
