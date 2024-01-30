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

function get_catalog(test::Bool)
    L  = [250.0]
    t  = [17.5]
    Lc = [15.0]

    set_L  = [400.0]
    set_t  = [40.0]
    set_Lc = [40.0]
    out = 0
    #loop here
    for i in [1]
        L = set_L[i]
        t = set_t[i]
        Lc = set_Lc[i]

        out = get_catalog(L,t,Lc, test = test)
    end

    return out
end


function get_catalog(L,t,Lc; test=true)
    if !test
        range_fc′ = 28.:1.:56.
        # range_as = 90.0:10.0:140  #[99.0, 140.0]
        range_as = [99.0,140.0]
        range_ec = 0.5:0.025:1.2
        range_fpe = (0.00:0.025:0.7) * 1860.0
    elseif test
        #test
        range_fc′ = 28.
        range_as = 140.0
        range_ec = 0.05
        range_fpe = 186.0
    else 
        println("Error Invalid test case")
        return nothing
    end

    nfc′ = length(range_fc′)
    nas  = length(range_as)
    nec  = length(range_ec)
    nfpe = length(range_fpe)
    ntotal = nfc′ * nas * nec * nfpe
    n = [nfc′, nas, nec, nfpe]


    results = Matrix{Float64}(undef, ntotal, 8 + 3) #L t Lc
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