using Dates
using DataFrames
using CSV
# using BenchmarkTools
using AsapSections
using Printf

include("Geometry/pixelgeo.jl")
# include("sectionproperties.jl")
include("Functions/postTensionedFunc.jl")
include("Functions/embodiedCarbon.jl")
include("Functions/capacities.jl")



# range_fc′ = 28.:7.:56.
# range_as = [99.0, 140.0]
# ec_max = 0.7
# range_ec = 0.5:0.1:ec_max
# range_fpe = (0.1:0.1:0.7) * 1860.0


#old version
# function map(idx_fc′, idx_as, idx_ec, idx_fpe;
#     range_fc′=range_fc′, range_as=range_as, range_ec=range_ec, range_fpe=range_fpe)
#     return idx_fpe +
#            length(range_fpe) * (idx_ec - 1) +
#            length(range_fpe) * length(range_ec) * (idx_as - 1) +
#            length(range_fpe) * length(range_ec) * length(range_as) * (idx_fc′ - 1)
# end

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

function get_catalog(test::Bool)
    L  = 250.0
    t  = 17.5
    Lc = 15.0
    return get_catalog(L,t,Lc, run_test = test)
end

function get_catalog(L,t,Lc; run_test=true)
    if !run_test
        range_fc′ = 28.:4.:56.
        range_as = 40.0:20.0:140#[99.0, 140.0]
        range_as = [99.0,140.0]
        range_ec = 0.05:0.025:1.2
        range_fpe = (0.00:0.025:0.7) * 1860.0
    else
        #test
        range_fc′ = 28.
        range_as = 99.0
        range_ec = 0.5
        range_fpe = 186.0
    end

    nfc′ = length(range_fc′)
    nas = length(range_as)
    nec = length(range_ec)
    nfpe = length(range_fpe)


    total_s = nfc′ * nas * nec * nfpe
    results = Matrix{Float64}(undef, total_s, 8)
    #we will loop through these three parameters and get the results.
    # with constant cross section properties.
    n = [nfc′, nas, nec, nfpe]
    for idx_fc′ in eachindex(range_fc′)
        for idx_as in eachindex(range_as)
            for idx_ec in eachindex(range_ec)
                for idx_fpe in eachindex(range_fpe)
                    fc′ = range_fc′[idx_fc′]
                    as = range_as[idx_as]
                    ec = range_ec[idx_ec]
                    fpe = range_fpe[idx_fpe]

                    pu, mu, vu, embodied = get_capacities(fc′, as, ec, fpe, L, t, Lc)
                    idx_all = [idx_fc′, idx_as, idx_ec, idx_fpe]

                    idx = mapping(n,idx_all)
                    results[idx,:] = [fc′, as, ec, fpe, pu, mu, vu, embodied]
                end
            end
        end
    end
    df = DataFrame(results , [ :fc′, :as,:ec,:fpe,:Pu,:Mu, :Vu, :carbon])
    df[!,:ID] = 1:total_s
    println("Got Catalog with $total_s outputs")
    return df# results # DataFrame(results)
end

#test
# results_test = get_catalog()
# results = get_catalog(100,10,10,run_test=false)
# # 11.147s , 575.72 MiB allocation
# date = Dates.today()
# time = Dates.now()

# CSV.write(joinpath(@__DIR__,"Outputs\\output_$date.csv"), results)

results = get_catalog(false)

CSV.write(joinpath(@__DIR__,"Catalogs/catalog_static.csv"), results)


# calcap(28., 99.0, 0.5, 1600.0)