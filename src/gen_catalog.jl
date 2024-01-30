using Dates
using DataFrames
using CSV
# using BenchmarkTools
using AsapSections
using Printf

include("Geometry/pixelgeo.jl")
include("Functions/embodiedCarbon.jl")
include("Functions/capacities.jl")

"""
Map an n dimentional vector into an index.
"""
function mapping(n::Vector{Int64}, idx::Vector{Int64})
    d = Vector{Int64}(undef, length(n))
    for i in eachindex(n) 
        if i == length(d)
            d[i] = mod(idx[i]+n[i]-1, n[i])+1
        else
            d[i] = (idx[i]-1)*prod(n[i+1:end])
        end
    end
    return sum(d)
end 


function get_catalog(test::Bool)::DataFrame
    pixel_sections = [ 205.0 35.0 30.0] #updated by Jenna Jan 2024
    out = DataFrame()
    for i in 1:size(pixel_sections)[1]
        L, t,Lc = pixel_sections[i,:]
        sub_catalog = get_catalog(L,t,Lc, test = test)
        out = vcat(out, sub_catalog)
    end
    return out
end

function get_catalog(L,t,Lc; test=true)::DataFrame
    if test
        #test
        range_fc′ = 28.
        range_as = 140.0
        range_ec = 0.05
        range_fpe = 186.0
    elseif !test
        
        range_fc′ = 28.:1.:56.
        range_as = [99.0*2,140.0*2] # x2 are for 2 ropes on 2 sides
        range_ec = 0.5:0.025:1.2
        range_fpe = (0.00:0.025:0.7) * 1860.0

    else 
        println("Error Invalid test case")
        return nothing
    end

    n = length.(range_fc′, range_as, range_ec, range_fpe)
    ntotal = prod(n)

    results = Matrix{Float64}(undef, prod(n), 8 + 3) #L t Lc
    #we will loop through these three parameters and get the results.
    # with constant cross section properties.
    for idx_fc′ in eachindex(range_fc′)
        for idx_as in eachindex(range_as)
            for idx_ec in eachindex(range_ec)
                for idx_fpe in eachindex(range_fpe)
                    fc′ = range_fc′[idx_fc′]
                    as = range_as[idx_as]
                    ec = range_ec[idx_ec]
                    fpe = range_fpe[idx_fpe]
                    #to do
                    #create a Section (ReinforcedConcreteSection or PixelFrameSection <: Section)
                    #then we do get_capacities(Section)
                    pu, mu, vu, embodied = get_capacities(fc′, as, ec, fpe, L, t, Lc)
                    idx_all = [idx_fc′, idx_as, idx_ec, idx_fpe]

                    idx = mapping(n,idx_all)
                    results[idx,:] = [fc′, as, ec, fpe, pu, mu, vu, embodied, L, t, Lc]
                end
            end
        end
    end
    df = DataFrame(results , [ :fc′, :as,:ec,:fpe,:Pu,:Mu, :Vu, :carbon, :L, :t,:Lc])
    df[!,:ID] = 1:ntotal
    println("Got Catalog with $ntotal outputs")
    return df# results # DataFrame(results)
end

#test
# results_test = get_catalog()
# results = get_catalog(100,10,10,test=false)
# # 11.147s , 575.72 MiB allocation
# date = Dates.today()
# time = Dates.now()

# CSV.write(joinpath(@__DIR__,"Outputs\\output_$date.csv"), results)

results = get_catalog(false)

CSV.write(joinpath(@__DIR__,"Catalogs/JAN23_5_catalog_static.csv"), results)


# calcap(28., 99.0, 0.5, 1600.0)