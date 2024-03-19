include("Functions//gen_catalog.jl")
# fc′_path = "src//Tables//fiber_with_extrapolation.csv"
fc′_path = "src//Tables//1503_fiber_for_test.csv"
results = get_catalog("default", fc′_path = fc′_path)
println("Saved Here")
CSV.write(joinpath(@__DIR__, "Catalogs/18_03_catalog_for_test.csv"), results)

#apply some design constraints to the created catalog.

#load some constraints
# hardware = CSV.read("Tables/Tension hardware list.csv", DataFrame)
