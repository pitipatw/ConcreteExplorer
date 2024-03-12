include("Functions//gen_catalog.jl")
fc′_path = "src//Tables//fiber_with_extrapolation.csv"
results = get_catalog("default", fc′_path = fc′_path)
println("Saved Here")
CSV.write(joinpath(@__DIR__, "Catalogs/MAR08_2_catalog_alltypes.csv"), results)
