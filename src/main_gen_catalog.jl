include("Functions//gen_catalog.jl")
fc′_path = "Tables//2302_fiber_with_extrapolation.csv"
results = get_catalog("default", fc′_path = fc′_path)
println("Here")
CSV.write(joinpath(@__DIR__, "Catalogs/FEB23_2_catalog_static.csv"), results)
