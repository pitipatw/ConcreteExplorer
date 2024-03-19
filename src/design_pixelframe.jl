using CSV, DataFrames, JSON
using Dates
# using ProgressBars
using UnPack
using Makie, GLMakie
using AsapToolkit #, kjlMakie

# set_theme!(kjl_light)
println(pwd())
include("Functions//definition.jl");
include("Functions//functions.jl");
# include("Functions/structuralelement.jl")
include("Functions//generalfunctions.jl");
include("Functions/preprocessing.jl");
include("Functions/pixel_design_opt.jl");
include("Functions//get_Deflection.jl");
include("Functions//interpolations.jl");
include("Functions/design_vis.jl")



"""
# PixelFrame Design Demo
"""

"""
###### To do list
1. Embodied carbon on the actual tendon profile
2. Turn everything into a gradient based optimization
"""

date = "18_03";
version = "01";
imagesavepath = "src//Images//Demo//" * date * "_";
println("Date: $date \nFile version $version")

# Load the design catalog
# Currently using version FEB6_4"""
# catalog = CSV.read("src/Catalogs/MAR08_1_catalog_alltypes.csv", DataFrame);
catalog = CSV.read("src//Catalogs//18_03_catalog_for_test.csv", DataFrame);
println(catalog[1:20, :]);

#load demands into a dictionary
# demand_path = joinpath(@__DIR__, "Demands/test_input_CISBAT_dataset.json");
first_part = [ "primary_beam", "secondary_beam","interior_column"] 
second_part = ["_occupancy_demands.json","_service_demands.json"] 
all_files = [i*j for i in first_part for j in second_part ] 

# demand_path = joinpath(@__DIR__, "16_03_test_demands.json");
# loaded_demands = JSON.parsefile(joinpath("src//test_demands.json"), dicttype =Dict{String, Vector{Union{Int64,Float64, String}}} );

demands = preprocessing(demand_path);

plot_distribution(demands, catalog)

# function designPixel(demands::DataFrame, catalog::DataFrame)

#Get feasible sections (subcatalog) for each demand point.
all_feasible_sections = filter_demands!(demands, catalog)

# #find keys that as 0 result 
# zero_result = []
# for k in keys(all_feasible_sections)
#     # println(k)
#     if length(all_feasible_sections[k]) == 1
#         if all_feasible_sections[k][1] == 0
#             push!(zero_result, k)
#         end
#     end
# end
# sort!(zero_result)

# #find elements of those sections 
# failed_elements = Dict{Int64,Vector{Float64}}()
# failed_type = []
# for k in zero_result
#     #k is the idx
#     e = demands[k, :e_idx]
#     if e ∉ failed_elements
#         push!(failed_elements, e)
#         push!(failed_type, demands[k, :type])
#     end
# end

# for i in zero_result[100:150]
#     println("###")
#     println(demands[i, :type])
#     println(demands[i, :load])
# end

# sort!(failed_elements)
# # sort!(failed_type, by = failed_elements)
# demands[130, :]
# for i in failed_elements
#     if demands[i*12, :type] == "columns"
#         println(i)
#         break
#     end
# end



# # Threads.nthreads()
# figures = Vector{Figure}(undef, size(demands)[1])
# #select a section to see the available designs
# # Threads.@threads for section_number in 1:size(demands)[1]
# for section_number in 1:size(demands)[1]
#     if section_number < 100
#         @show section_number
#     end
#     # section_number = 1

#     figure_check_section = Figure(size=(500, 500))
#     ax_1 = Axis(figure_check_section[1, 1], xlabel="Moment [kNm]", ylabel="Shear [kN]", title=string(section_number))
#     ax_2 = Axis(figure_check_section[2, 1], xlabel="dps", ylabel="T")

#     scatter!(ax_1, catalog[all_feasible_sections[section_number], :Mu], catalog[all_feasible_sections[section_number], :Vu], color=catalog[all_feasible_sections[section_number], :fc′], colorrange=extrema(catalog[!, :fc′]))
#     scatter!(ax_2, catalog[all_feasible_sections[section_number], :dps], catalog[all_feasible_sections[section_number], :T], color=catalog[all_feasible_sections[section_number], :fc′], colorrange=extrema(catalog[!, :fc′]))

#     scatter!(ax_1, demands[section_number, :mu], demands[section_number, :vu], marker='x')
#     type_map = Dict("primary" => 3, "secondary" => 3, "columns" => 2)
#     scatter!(ax_2, demands[section_number, :ec_max], getindex.(Ref(type_map), demands[section_number, :type]), marker='x')

#     figures[section_number] = figure_check_section
#     # figure_check_section
#     # @show imagesavepath * "figure_check_section" * string(section_number) * ".png"
#     # save(imagesavepath * "figure_check_section" * string(section_number) * ".png", figure_check_section)
# end

# for i in eachindex(figures)
#     save(imagesavepath * "figure_check_section" * string(i) * ".png", figures[i])
# end

# for i in 1:length(all_feasible_sections)
#     print("Section $i")
#     println(length(all_feasible_sections[i]))
# end

elements_designs, elements_to_sections, sections_to_designs, skipped_elements = find_optimum(all_feasible_sections, demands);
# println(elements_designs)

max_primary = [1, 0]
max_secondary = [1, 0]
max_column = [1, 0]

for i in eachindex(elements_designs)
    # @show i
    if i ∉ skipped_elements
        #get a section (the load will be the same)
        s = elements_to_sections[i][1]
        
        type = demands[s, :type]
        load = demands[s, :load]
        element = demands[s, :e_idx]
        if type == "primary"
            if max_primary[2] < load
                @show max_primary = [element, load]
            end
        elseif type == "secondary"
            if max_secondary[2] < load
                @show max_secondary = [element, load]
            end
        elseif type == "columns"
            if max_column[2] < load
                @show max_column = [element, load]
            end
        end
    end
end

@show max_primary
@show max_secondary
@show max_column




# elements_designs_fielded = Vector{Dict{String,Real}}()
# # for i in eachindex(elements_designs)
# open("src//Results//designs_results_15_03.json", "w") do f
#     JSON.print(f, elements_designs)
# end


# open("src//Results//sections_to_designs_15_03.json", "w") do f
#     JSON.print(f, sections_to_designs)
# end

using Makie, GLMakie, CairoMakie
using JSON
using DataFrames, CSV

# designs = JSON.parsefile(joinpath(@__DIR__, "Results//designs_results_08_03.json"), dicttype=Dict{String,Vector{Vector{Float64}}});
# mapping_strings = ["ID", "fc′", "dosage", "fR1", "fR3", "as", "dps", "fpe","fps", "Pu", "Mu", "Vu,", "carbon", "L", "t", "Lc", "T", "catalog_id", "max_dps", "min_dps"]
# for i in eachindex(elements_designs)
#     for s in eachindex(elements_designs[i])
#         global_s_index = elements_to_sections[i][s]
#         println(global_s_index)
#         # @show  Dict(mapping_strings .=> elements_designs[parse(Int64,i)][s])
#         # @show elements_designs_fielded[global_s_index]
#         push!(elements_designs_fielded, Dict(mapping_strings .=> vcat(global_s_index - 1, elements_designs[i][s])))
#         # elements_designs_fielded[global_s_index] =
#     end
# end

# open("src/Results/12_03_designs_results_fielded.json", "w") do f
#     JSON.print(f, elements_designs_fielded)
# end

# ne = length(designs)
# println("There are $ne elements.")




#test 
plot_element(max_primary[1], elements_designs)
plot_element(max_secondary[1], elements_designs)
plot_element(max_column[1], elements_designs)



for e_idx in [max_primary[1], max_secondary[1], max_column[1]]
    f = plot_element(e_idx, elements_designs)
    save("src/Results/"*date * "test_$e_idx.png", f)
end

primary_design   = elements_designs[max_primary[1]]
secondary_design = elements_designs[max_secondary[1]]
column_design    = elements_designs[max_column[1]]


mymapping = keys(primary_design[1])
output = DataFrame([name => [] for name in mymapping])
for i in 1:10
    push!(output,primary_design[i])    
end
for i in 1:11
    push!(output,secondary_design[i])    
end
for i in 1:12
    push!(output,column_design[i])    
end

open("src//Results//16_03_design_results_fortest.json", "w") do f
    JSON.print(f, output)
end

@show 1
# for i in 1:19
#     f = plot_element(i, designs)
#     save("src/Results/Results_" * date * "/$i.png", f)
# end

# #summarize the result. 
# function get_design_properties(sections_to_designs::Dict{Int64,Vector{Float64}}, idx::Int64)
#     output = Vector{Float64}(undef, length(sections_to_designs))
#     for i in eachindex(sections_to_designs)
#         @show i
#         output[i] = sections_to_designs[i][idx]
#     end
#     return output
# end



# all_fc′ = get_design_properties(sections_to_designs, 1)
# all_dosage = get_design_properties(sections_to_designs, 2)
# all_fR1 = get_design_properties(sections_to_designs, 3)
# all_fR3 = get_design_properties(sections_to_designs, 4)
# all_as = get_design_properties(sections_to_designs, 5)
# all_dps = get_design_properties(sections_to_designs, 6)
# all_fpe = get_design_properties(sections_to_designs, 7)
# stack_name = hcat(string.(all_fc′, "_", all_fR1, "_", all_fR3, "_", all_dosage))
# @show unique_stack_name = unique(stack_name)
# MK_file_prep = DataFrame(:fc′ => all_fc′, :fR1 => all_fR1, :fR3 => all_fR3, :dosage => all_dosage)

# csv_fc′ = Vector{Float64}()
# csv_dosage = Vector{Real}()
# csv_fR1 = Vector{Float64}()
# csv_fR3 = Vector{Float64}()

# for i in 1:length(unique_stack_name)
#     @show vals = parse.(Float64, split(unique_stack_name[i], "_"))
#     @show typeof(vals[1])
#     @show typeof(vals[2])
#     @show typeof(vals[3])
#     push!(csv_fc′, vals[1])
#     push!(csv_fR1, vals[2])
#     push!(csv_fR3, vals[3])
#     push!(csv_dosage, vals[4])

# end

# csv_output = DataFrame(:fc′ => csv_fc′, :dosage => csv_dosage, :fR1 => csv_fR1, :fR3 => csv_fR3)
# CSV.write(joinpath(@__DIR__, date * "_mix_specs.csv"), csv_output)



# #each pair, plots them dots and x and a line connecting them together. 

# f_final = Figure(size=(500, 500))
# ax1 = Axis(f_final[1, 1], xlabel="Moment [kNm]", ylabel="Shear [kN]", title="Demands vs Designs")
# demand_points = hcat(demands[!, :mu], demands[!, :vu])
# design_points = hcat(get_design_properties(sections_to_designs, 9), get_design_properties(sections_to_designs, 10))
# for i in 1:size(demand_points)[1]
#     x1 = demand_points[i, 1]
#     y1 = demand_points[i, 2]
#     x2 = design_points[i, 1]
#     y2 = design_points[i, 2]
#     @assert x2 > x1
#     @assert y2 > y1 i
#     u = x2 - x1
#     v = y2 - y1
#     arrows!([x1], [y1], [u], [v], arrowsize=5)
# end
# scatter!(ax1, demand_points[:, 1], demand_points[:, 2], color=:red, markersize=10)
# scatter!(ax1, design_points[:, 1], design_points[:, 2], color=all_fc′, market_size=10)

# f_final
