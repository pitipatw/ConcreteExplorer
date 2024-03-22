using CSV, DataFrames, JSON
using Dates
# using ProgressBars
using UnPack
using Makie, GLMakie
using AsapToolkit #, kjlMakie

# set_theme!(kjl_light)
println(pwd())


include("Functions/catalog_vis.jl")
# include("Functions/structuralelement.jl")
include("Functions/preprocessing.jl");
include("Functions/pixel_design_search.jl");
include("Functions/design_vis.jl")

include("Functions//Deflection//get_Deflection.jl");
include("Functions//Deflection//interpolations.jl");
include("Functions//Deflection//definition.jl");
include("Functions//Deflection//functions.jl");
include("Functions//Deflection//generalfunctions.jl");



"""
# PixelFrame Design Demo
"""

"""
###### To do list
1. Embodied carbon on the actual tendon profile
2. Turn everything into a gradient based optimization
"""


date = "21_03";
version = "01";
date *= "_"*version
save_directory = "src//Results//"*date*"//";
save_image_directory = save_directory*"Images//";


if !ispath(save_directory)
    mkdir(save_directory)
end
additional_path = "20_03_FullScaleTestDemands//" #for more specific task
save_directory *=additional_path
if !ispath(save_directory)
    mkdir(save_directory)
end

if !ispath(save_image_directory)
    mkdir(save_image_directory)
end

if !ispath(save_directory*additional_path)
    mkdir(save_directory*additional_path)
end

println("Current date: $date.\nFile version $version\nFiles will be saved at \n :$save_directory  \n :$save_image_directory")

# Load the design catalog
# Currently using version FEB6_4"""
# catalog = CSV.read("src/Catalogs/MAR08_1_catalog_alltypes.csv", DataFrame);
catalog_path = "src//Catalogs//18_03_catalog_for_test.csv"
catalog = CSV.read(catalog_path, DataFrame);
println("Catalog was loaded from $catalog_path")
println(catalog[1:5, :]);
f_catalog = catalog_vis(catalog_path)
#load demands into a dictionary
# demand_path = joinpath(@__DIR__, "Demands/test_input_CISBAT_dataset.json");
first_part = [ "8m_primary_beam", "8m_secondary_beam","8m_interior_column"] ;
second_part = ["_occupancy_demands"];#,"_service_demands"] 



all_files = [i*j for i in first_part for j in second_part ] ;

# for i in eachindex(all_files) 
    i = 1 
    filename = all_files[i]
    demand_path = joinpath(@__DIR__, "Demands//"*additional_path*filename*".json");

    demands = preprocessing(demand_path, save_image_directory);

    plot_distribution(demands, catalog, save_image_directory);

    all_feasible_sections = filter_demands!(demands, catalog)

    plot_feasible_sections(all_feasible_sections, save_image_directory)


    elements_designs, elements_to_sections, sections_to_designs, skipped_elements = search_design(all_feasible_sections, demands);

    for i in keys(elements_designs) 
        f = plot_element(i, elements_designs, elements_to_sections)
        save(save_image_directory*"_output_from_"*filename*"_$i.png", f)
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
# end