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



"""
# PixelFrame Design Demo
"""

"""
###### To do list
1. Embodied carbon on the actual tendon profile
2. Turn everything into a gradient based optimization
"""

date = "13_03";
version = "01";
imagesavepath = "src//Images//Demo//" * date * "_";
println("Date: $date \nFile version $version")

# Load the design catalog
# Currently using version FEB6_4"""
catalog = CSV.read("src/Catalogs/MAR08_1_catalog_alltypes.csv", DataFrame);
println(catalog[1:20, :]);

#load demands into a dictionary
# demand_path = joinpath(@__DIR__, "Demands/test_input_CISBAT_dataset.json");
demand_path = joinpath(@__DIR__, "Demands/0308_test_building_2_occupancy.json");

demands = preprocessing(demand_path);

plot_distribution(demands, catalog)

#Get feasible sections (subcatalog) for each demand point.
all_feasible_sections = filter_demands!(demands, catalog)

#select a section to see the available designs
for section_number in 1:size(demands)[1]
    # section_number = 1
    figure_check_section = Figure(size=(500, 500))
    ax_1 = Axis(figure_check_section[1, 1], xlabel="Moment [kNm]", ylabel="Shear [kN]", title=string(section_number))
    ax_2 = Axis(figure_check_section[2, 1], xlabel="dps", ylabel="T")

    for current_mid_design in all_feasible_sections[section_number]
        scatter!(ax_1, catalog[current_mid_design, :Mu], catalog[current_mid_design, :Vu], color=catalog[current_mid_design, :fc′], colorrange=extrema(catalog[!, :fc′]))
        scatter!(ax_2, catalog[current_mid_design, :dps], catalog[current_mid_design, :T], color=catalog[current_mid_design, :fc′], colorrange=extrema(catalog[!, :fc′]))
    end
    scatter!(ax_1, demands[section_number, :mu], demands[section_number, :vu], marker='x')
    type_map = Dict("primary" => 3, "secondary" => 3, "columns" => 2)
    scatter!(ax_2, demands[section_number, :ec_max], getindex.(Ref(type_map), demands[section_number, :type]), marker='x')

    figure_check_section
    @show imagesavepath * "figure_check_section" * string(section_number) * ".png"
    save(imagesavepath * "figure_check_section" * string(section_number) * ".png", figure_check_section)
end

for i in 1:length(all_feasible_sections)
    println("Section $i")
    println(length(all_feasible_sections[i]))
end

elements_designs, elements_to_sections, sections_to_designs = find_optimum(all_feasible_sections, demands)
println(elements_designs)


elements_designs_fielded = Vector{Dict{String,Real}}()
# for i in eachindex(elements_designs)
open("src//Results//designs_results_08_03.json", "w") do f
    JSON.print(f, elements_designs)
end


open("src//Results//sections_to_designs_08_03.json", "w") do f
    JSON.print(f, sections_to_designs)
end

using Makie, GLMakie, CairoMakie
using JSON
using DataFrames, CSV

designs = JSON.parsefile(joinpath(@__DIR__, "Results//designs_results_08_03.json"), dicttype=Dict{String,Vector{Vector{Float64}}});
mapping_strings = ["ID", "fc′", "dosage", "fR1", "fR3", "as", "dps", "fpe", "Pu", "Mu", "Vu,", "carbon", "L", "t", "Lc", "T", "catalog_id", "max_dps", "min_dps"]
for i in eachindex(elements_designs)
    for s in eachindex(elements_designs[i])
        global_s_index = elements_to_sections[i][s]
        println(global_s_index)
        # @show  Dict(mapping_strings .=> elements_designs[parse(Int64,i)][s])
        # @show elements_designs_fielded[global_s_index]
        push!(elements_designs_fielded, Dict(mapping_strings .=> vcat(global_s_index - 1, elements_designs[i][s])))
        # elements_designs_fielded[global_s_index] =
    end
end

open("src/Results/12_03_designs_results_fielded.json", "w") do f
    JSON.print(f, elements_designs_fielded)
end

ne = length(designs)
println("There are $ne elements.")

function plot_element(element_number::Int64, designs::Dict;
    L::Float64=250.0)

    sections = elements_to_sections[element_number]
    element_number = string(element_number)

    L = designs[element_number][1][12]
    @show set_fc′ = [i[1] for i in designs[element_number]]
    tendon_profile = [i[6] for i in designs[element_number]]
    axial_capacity = [i[8] for i in designs[element_number]]
    moment_capacity = [i[9] for i in designs[element_number]]
    shear_capacity = [i[10] for i in designs[element_number]]


    axial_demand = [demands[i, "pu"] for i in sections]
    moment_demand = [demands[i, "mu"] for i in sections]
    shear_demand = [demands[i, "vu"] for i in sections]


    #plot center around x = 0 
    @show n = length(tendon_profile)
    @show length(axial_demand)
    @show res = mod(n, 2) * 250
    xmax = div(n, 2) * 500
    x_1 = -xmax + 250 - res
    x_n = xmax - 250 - res
    # @show x_range = (-xmax+res:500:xmax-res)
    @show x_range = x_1:500:x_n

    #create a bands (polygon of possible tendon profile)
    tendon_points = Matrix{Int64}(undef, 2, 2 * n)
    for i in 1:n
        tendon_points[:, i] = [500 * (i - 1) - xmax, -designs[element_number][i][17]]
    end
    for i in 1:n
        tendon_points[:, 2*n-i+1] = [500 * (i - 1) - xmax, -designs[element_number][i][18]]
    end

    f1 = Figure(size=(1200, 600))
    g = f1[1, 1] = GridLayout()
    axs_design = Axis(g[1, 1], title="Element $element_number", titlesize=20,
        aspect=DataAspect(),
        limits=(x_range[1] - 100, x_range[end] + 600, -2 * L, 0.8 * L),
        yticks=div(L, 100)*100:-100:-div(L, 100)*105,
        # yminorticks = IntervalsBetween(2),
        # yminorgridvisible = true,
        ylabel="y"
    )

    for i in 1:n
        points = [-xmax + res + (i - 1) * 500, -L, 500, L * 1.5]
        poly!(axs_design, Rect(points...), color=set_fc′[i], colorrange=(20, 100), colormap=:grays)
    end
    Colorbar(f1[1, 2], ticks=[20, 30, 40, 50, 60, 70, 80, 90, 100], colorrange=(40, 80), colormap=cgrad(:grays, 3, categorical=true, rev=true))
    #plot section based on concrete strength 
    # concrete = poly!(axs_design,Rect( [-xmax+res, -L, (n-1)*500, L*1.5]...), color = (:grey,0.2))
    # concrete = poly!(axs_design,Rect( -xmax+res, -L, (n-1)*500, L*1.5), color = (:grey,0.2))

    tendon = poly!(axs_design, tendon_points, color=:skyblue, alpha=0.1, transparent=true)
    tendon_pts = scatter!(axs_design, tendon_points)
    tendon = lines!(axs_design, x_range, -tendon_profile)

    axs_axial = Axis(g[2, 1], aspect=10,
        limits=(x_range[1], x_range[end] + 600, nothing, nothing), ylabel="Axial [kN]",
    )

    axs_moment = Axis(g[3, 1], aspect=10,
        limits=(x_range[1], x_range[end] + 600, nothing, nothing), ylabel="Moment [kNm]",
    )
    axs_shear = Axis(g[4, 1], aspect=10,
        limits=(x_range[1], x_range[end] + 600, nothing, nothing), ylabel="Shear [kN]",
    )
    hidexdecorations!(axs_design, grid=false)
    hidexdecorations!(axs_axial, grid=false)
    hidexdecorations!(axs_moment, grid=false)

    lines!(axs_axial, x_range, axial_capacity, color=:red)
    lines!(axs_moment, x_range, moment_capacity, color=:blue)
    lines!(axs_shear, x_range, shear_capacity, color=:green)

    lines!(axs_axial, x_range, axial_demand, linestyle=:dash, color=:red)
    lines!(axs_moment, x_range, moment_demand, linestyle=:dash, color=:blue)
    lines!(axs_shear, x_range, shear_demand, linestyle=:dash, color=:green)


    # for (i, label) in enumerate(["Axial [kN]", "Moment [kNm]", "Shear [kN]"])
    #     Box(g[i, 2], color = :gray90)
    #     Label(g[i,2], label, rotation = pi/2, tellheight = false)
    # end

    rowgap!(g, 10)

    yspace = maximum(tight_yticklabel_spacing!, [axs_axial, axs_shear, axs_moment]) + 10

    axs_axial.yticklabelspace = yspace
    axs_moment.yticklabelspace = yspace
    axs_shear.yticklabelspace = yspace

    return f1
end


#test 
plot_element(1, designs)

for i in 1:19
    f = plot_element(i, designs)
    save("src/Results/Results_" * date * "/$i.png", f)
end

#summarize the result. 
function get_design_properties(sections_to_designs::Dict{Int64,Vector{Float64}}, idx::Int64)
    output = Vector{Float64}(undef, length(sections_to_designs))
    for i in eachindex(sections_to_designs)
        @show i
        output[i] = sections_to_designs[i][idx]
    end
    return output
end



all_fc′ = get_design_properties(sections_to_designs, 1)
all_dosage = get_design_properties(sections_to_designs, 2)
all_fR1 = get_design_properties(sections_to_designs, 3)
all_fR3 = get_design_properties(sections_to_designs, 4)
all_as = get_design_properties(sections_to_designs, 5)
all_dps = get_design_properties(sections_to_designs, 6)
all_fpe = get_design_properties(sections_to_designs, 7)
stack_name = hcat(string.(all_fc′, "_", all_fR1, "_", all_fR3, "_", all_dosage))
@show unique_stack_name = unique(stack_name)
MK_file_prep = DataFrame(:fc′ => all_fc′, :fR1 => all_fR1, :fR3 => all_fR3, :dosage => all_dosage)

csv_fc′ = Vector{Float64}()
csv_dosage = Vector{Real}()
csv_fR1 = Vector{Float64}()
csv_fR3 = Vector{Float64}()

for i in 1:length(unique_stack_name)
    @show vals = parse.(Float64, split(unique_stack_name[i], "_"))
    @show typeof(vals[1])
    @show typeof(vals[2])
    @show typeof(vals[3])
    push!(csv_fc′, vals[1])
    push!(csv_fR1, vals[2])
    push!(csv_fR3, vals[3])
    push!(csv_dosage, vals[4])

end

csv_output = DataFrame(:fc′ => csv_fc′, :dosage => csv_dosage, :fR1 => csv_fR1, :fR3 => csv_fR3)
CSV.write(joinpath(@__DIR__, date * "_mix_specs.csv"), csv_output)



#each pair, plots them dots and x and a line connecting them together. 

f_final = Figure(size=(500, 500))
ax1 = Axis(f_final[1, 1], xlabel="Moment [kNm]", ylabel="Shear [kN]", title="Demands vs Designs")
demand_points = hcat(demands[!, :mu], demands[!, :vu])
design_points = hcat(get_design_properties(sections_to_designs, 9), get_design_properties(sections_to_designs, 10))
for i in 1:size(demand_points)[1]
    x1 = demand_points[i, 1]
    y1 = demand_points[i, 2]
    x2 = design_points[i, 1]
    y2 = design_points[i, 2]
    @assert x2 > x1
    @assert y2 > y1 i
    u = x2 - x1
    v = y2 - y1
    arrows!([x1], [y1], [u], [v], arrowsize=5)
end
scatter!(ax1, demand_points[:, 1], demand_points[:, 2], color=:red, markersize=10)
scatter!(ax1, design_points[:, 1], design_points[:, 2], color=all_fc′, market_size=10)

f_final
