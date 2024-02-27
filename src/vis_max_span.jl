using JSON, DataFrames
using Makie, GLMakie

output = JSON.parsefile(joinpath(@__DIR__,"Results//26_02_span.json"), dicttype = Dict{Any,Any}) #Dict{String, Vector{Dict{String, Union{Real, String}}}});

#visualziaing the result. 

#count output
n_output = 0
for i in 1:100
    for configs in output[string(i)] 
        n_output +=1 
    end
end

n_output


f_output = Figure(size = (500,500))
ax1 = Axis(f_output[1,1], title = "ID vs Sections")
ax2 = Axis(f_output[1,2], title = "Span vs Load")
ax3 = Axis(f_output[1,3], title = "L vs Load")
ax6 = Axis(f_output[3,3], title = "L vs Span")


loads = Vector{Float64}(undef, n_output)
idx = Vector{Float64}(undef, n_output)
span = Vector{Int64}(undef, n_output)
set_L = Vector{Float64}(undef, n_output)
set_Lc = Vector{Float64}(undef, n_output)
set_t = Vector{Float64}(undef, n_output)

check = Vector{Symbol}(undef, n_output)
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


        stat  = configs["control"]
        if stat == "Shear control"
            check[count] = :red
        elseif stat == "Moment control"
            check[count] = :blue 
        end
    end
end

scatter!(ax1, idx, loads, color = check)
scatter!(ax2, span, loads, color = check)
scatter!(ax3, set_L, loads, color = check)
scatter!(ax6, set_L, span, color = check)



f_output


f_all_sections

# find capacity of the section 



