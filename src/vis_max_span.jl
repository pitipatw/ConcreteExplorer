using JSON, DataFrames
using Makie, GLMakie

output = JSON.parsefile(joinpath(@__DIR__,"Results//27_02_span.json"), dicttype = Dict{Any,Any}) #Dict{String, Vector{Dict{String, Union{Real, String}}}});

#visualziaing the result. 

#count output
global n_output = 0
for i in 1:100
    for configs in output[string(i)] 
        global n_output +=1 
    end
end

n_output


f_output = Figure(size = (1500,1500))
ax1 = Axis(f_output[1,1], title = "ID vs Loads [kN]", xlabel = "ID", ylabel = "Load [kN]")
ax2 = Axis(f_output[1,2], title = "Span vs Load", xlabel = "Span [mm]", ylabel = "Load [kN]")
ax3 = Axis(f_output[1,3], title = "L vs Load", xlabel = "L [mm]", ylabel = "Load [kN]")
ax4 = Axis(f_output[2,1], title = "L vs t on Load", xlabel = "L [mm]", ylabel = "t [mm]")
ax5 = Axis3(f_output[2,2], title = "L vs Span vs Load", xlabel = "L [mm]", ylabel = "Span [mm]", zlabel = "Load [kN]")
ax6 = Axis(f_output[2,3], title = "L vs dps", xlabel = "L [mm]", ylabel = "dps [mm]")

ax9 = Axis3(f_output[3,3], title = "L vs Load vs span (color:as)", xlabel = "L [mm]", ylabel = "Load [kN]", zlabel = "span [mm]")

loads = Vector{Float64}(undef, n_output)
idx = Vector{Float64}(undef, n_output)
span = Vector{Int64}(undef, n_output)
set_L = Vector{Float64}(undef, n_output)
set_Lc = Vector{Float64}(undef, n_output)
set_t = Vector{Float64}(undef, n_output)
set_dps = Vector{Float64}(undef, n_output)
set_as = Vector{Float64}(undef, n_output)
color = Vector{Union{Symbol, Tuple{Symbol,Float64}}}(undef, n_output)
marker = Vector{Symbol}(undef, n_output)
global count = 0
for i in 1:100
    # L,t,Lc = all_sections[i,:]
    for configs in output[string(i)]
        global count +=1 
        idx[count] = i
        loads[count] = configs["load"]
        span[count] = configs["span"]
        set_L[count] = configs["L"]
        set_Lc[count] = configs["Lc"]
        set_t[count] = configs["t"]
        set_dps[count] = configs["dps"]
        set_as[count] = configs["as"]

        stat  = configs["control"]
        if stat == "Shear control"
            color[count] = (:grey, 0.1)
            marker[count] = :cross
        elseif stat == "Moment control"
            color[count] = :blue 
            marker[count] = :circle
        end
    end
end

#filter moment only

filter_criteria = color .== :blue .&& set_dps .== set_L


scatter!(ax1, idx, loads, color = color)
scatter!(ax2, span, loads, color = color)
scatter!(ax3, set_L, loads, color = color)
# scatter!(ax4, set_L, set_t, color = loads, marker = marker)
scatter!(ax4, set_L[filter_criteria], set_t[filter_criteria], color = loads[filter_criteria], marker = marker[filter_criteria])
scatter!(ax5, set_L[filter_criteria], span[filter_criteria], loads[filter_criteria],color = loads[filter_criteria], marker = marker[filter_criteria])
scatter!(ax6, set_L[filter_criteria], set_dps[filter_criteria], color = loads[filter_criteria])
# lines!(ax6, set_L[filter_criteria], set_L[filter_criteria], color = :orange)
s9 = scatter!(ax9, set_L[filter_criteria],loads[filter_criteria], span[filter_criteria], color = set_as[filter_criteria])
Colorbar(f_output[4,3],s9, vertical = false, tellwidth = false)


f_output


#filtered output 
# filtered_output = DataFrame(output)



|