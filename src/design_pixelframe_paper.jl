using CSV, DataFrames, JSON
using Dates
# using ProgressBars
using UnPack
using Makie, GLMakie, CairoMakie
using AsapToolkit #, kjlMakie

# set_theme!(kjl_light)
println(pwd())
# include("Functions//definition.jl");
# include("Functions//functions.jl");
# include("Functions/structuralelement.jl")
# include("Functions//generalfunctions.jl");
include("Functions/preprocessing.jl");
include("Functions/pixel_design_search_paper.jl");
# include("Functions//get_Deflection.jl");
# include("Functions//interpolations.jl");
include("Functions/design_vis.jl");



"""
# PixelFrame Design Demo
"""

"""
###### To do list
1. Embodied carbon on the actual tendon profile
2. Turn everything into a gradient based optimization
"""

date = "27_03";
version = "test_";
version = "building_test"
date = date*"_"*version

println("Date: $date \nFile version $version")

# Load the design catalog
# Currently using version FEB6_4"""
# catalog = CSV.read("src/Catalogs/MAR08_1_catalog_alltypes.csv", DataFrame);
catalog = CSV.read("src//Catalogs//20_03_catalog_for_paper.csv", DataFrame);
# catalog = CSV.read("src//Catalogs//18_03_catalog_for_test.csv", DataFrame);
println(catalog[1:20, :]);
catalog[!,:ID] = 1:size(catalog)[1]

#load demands into a dictionary
# demand_path = joinpath(@__DIR__, "Demands/test_input_CISBAT_dataset.json");
first_part = [ "primary_beam", "secondary_beam","interior_column"] 
second_part = ["_occupancy_demands"]#,"_service_demands"] 

folder_path = "FullScaleTestDemands//"
all_files = [i*j for i in first_part for j in second_part ] 

folder_path = ""
all_files = ["0902_three story demands_service"]
for i in eachindex(all_files) 
    # i = 1 
    filename = all_files[i]
    @show demand_path = joinpath(@__DIR__, "Demands//"* folder_path*filename*".json");
    imagesavepath = "src//Images//Demo//" * date * "_$filename"*"_";

    # demand_path = joinpath(@__DIR__, "16_03_test_demands.json");
    # loaded_demands = JSON.parsefile(joinpath("src//test_demands.json"), dicttype =Dict{String, Vector{Union{Int64,Float64, String}}} );

    demands = preprocessing(demand_path, imagesavepath);

    plot_distribution(demands, catalog, imagesavepath)

    # function designPixel(demands::DataFrame, catalog::DataFrame)

    #Get feasible sections (subcatalog) for each demand point.
    all_feasible_sections = filter_demands!(demands, catalog)

    for c in [1,2,3,4,5]
    end 

    elements_designs, elements_to_sections, sections_to_designs, skipped_elements = search_design(all_feasible_sections, demands);

    if !(ispath("src//Results//"*date))
        mkpath(date)
    end

    for i in keys(elements_designs) 
        f = plot_element(i, elements_designs, elements_to_sections)
        save(save_image_directory*date*"_"*filename* "$i.png", f)
    end 


    mymapping = keys(elements_designs[1][1])
    output = DataFrame([name => [] for name in mymapping])

    for i in keys(elements_designs)
        for j in sort!(collect(keys(elements_designs[i])))
            @show elements_designs[i][j]
        push!(output,elements_designs[i][j])    
        end
    end

    open(save_directory*"_output_from_"*filename*".json", "w") do f
        JSON.print(f, output)
    end
end