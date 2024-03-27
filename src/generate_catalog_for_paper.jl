include("Functions//gen_catalog.jl")
include("Functions//catalog_vis.jl")
# fc′_path = "src//Tables//fiber_with_extrapolation.csv"
# fc′_path = "src//Tables//18_03_fc_fiber_for_paper.csv"

fc′_path = "src//Tables//18_03_fc_fiber_for_paper.csv"
results = get_catalog("myResults", fc′_path = fc′_path)
# results = get_catalog("myResults", fc′_path = fc′_path)
println("Saved Here")
catalog_path = "Catalogs/27_03_catalog_for_paper.csv"
CSV.write(joinpath(@__DIR__,catalog_path), results)
f_catalog = catalog_vis(joinpath(@__DIR__,catalog_path))
f_catalog
#apply some design constraints to the created catalog.

#load some constraints
# hardware = CSV.read("Tables/Tension hardware list.csv", DataFrame)
