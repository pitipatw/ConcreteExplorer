# module PostTen
using JSON, HTTP
using CSV, DataFrames
using Dates


# include("pixelgeo.jl") #generating Pixel geometries
# include("sectionproperties.jl")
# include("calstr.jl") #calculating strength
# include("getterrain.jl")
# include("connection.jl")
# include("connection_dummy.jl")
# include("plotjson.jl")


"""
cin format
fc', as, ec, fpe, pu, mu, vu, embodied
"""
catalog = get_catalog()
catalog = DataFrame(x = 1:100,y = 1:2:200)
N = lengtg(catalog)
demands = JSON.parsefile("src/ReinforcedConcrete/JSON/demands_dummy.json")
# demands = (CSV.read("JSON/demands_dummy.csv", DataFrame))
designs = Dict{Int64,Vector{Int64}}() #mapping demands index into design index in the category
for d_idx in eachindex(demands) 
    pu = demands[d_idx]["pu"]
    mu = demands[d_idx]["mu"]
    # vu = demands[d_idx]["vu"] # neglect shear

    design_i = Vector{Int64}()
    filter1 = catalog[!,"pu"] .== pu
    filter2 = catalog[!,"mu"] .== mu

    design_i = catalog[ filter1 && filter2, :]

    designs[d_idx] = copy(design_i)
end

#box plot to visualize the design for each section.

#select the best design

# sort each design by gwp from low to high, and get the first one (lowest)\
# design
# gwpi = catalog[ design[i to n]] 
# sort by gwp.

# select gwpi[1]


# output designs_opt = Dict{Int64,Int64}
# have to work with the catalog.
# #

# visualize.
# plot by section?  needs to know section information 
# what section number maps to what elements. 
# We can do 1 section for the whole beam. 

# ###########

# fc′_range = [28.0 , 35.0] # or [28.0, 28.0] if you only want 28 MPa.
# # this should show how having a choie is better or not.
# loop designs
# for each one, get fc′ in the range. 

#     still select the least ec. 


