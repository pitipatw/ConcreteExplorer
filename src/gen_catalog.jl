using Dates
using DataFrames
using CSV
# using BenchmarkTools
using AsapSections
using Printf

include("Functions/Geometry/pixelgeo.jl")
include("Functions/embodiedCarbon.jl")
include("Functions/capacities.jl")

"""
Map an n dimentional vector into an index.
"""
function mapping(n::Vector{Int64}, idx::Vector{Int64})::Int64
    d = Vector{Int64}(undef, length(n))
    for i in eachindex(n)
        if i == length(d)
            d[i] = mod(idx[i] + n[i] - 1, n[i]) + 1
        else
            d[i] = (idx[i] - 1) * prod(n[i+1:end])
        end
    end
    return sum(d)
end


function get_catalog(test::Bool)::DataFrame
    pixel_sections = [205.0 35.0 30.0] #updated by Jenna Jan 2024
    out = DataFrame()
    for i in 1:size(pixel_sections)[1]
        L, t, Lc = pixel_sections[i, :]
        sub_catalog = get_catalog(L, t, Lc, test=test)
        out = vcat(out, sub_catalog)
    end
    return out
end

function get_catalog(L, t, Lc; test=true)::DataFrame
    if test #depreciated
        #test
        range_fc′ = 28.0
        range_as = 140.0
        range_dps = 50.0
        range_fpe = 186.0
    elseif !test
        fc_fiber = CSV.read("src//fc_fiber.csv", DataFrame)
        #These come in pair
        range_fc′= convert(Array{Float64}, fc_fiber[!, :strength])
        range_fR1 = fc_fiber[!, :fR1]
        range_fR3 = fc_fiber[!, :fR3]
        @assert length(range_fc′) == length(range_fR1) == length(range_fR3)

        range_as = [99.0 * 2, 140.0 * 2] # x2 are for 2 ropes on 2 sides
        range_dps = vcat(0.0, 50:10:300) #mm either in the middle, or 5cm shifted.
        range_fpe = (0.00:0.025:0.7) * 1860.0 #MPa
        range_type = [3.0, 2.0, 4.0] #PixelFrame configuration -> Y = 3 ,X2 = 2, X4 = 4.
    else
        println("Error Invalid test case")
        return nothing
    end

    n = length.([range_fc′, range_as, range_dps, range_fpe, range_type])
    ntotal = prod(n)
    #Pre allocating results
    results = Matrix{Float64}(undef, prod(n), 9 + 2 + 3) #L t Lc
    #we will loop through these three parameters and get the results.
    # with constant cross section properties.
    for idx_type in eachindex(range_type)
        for idx_fc′ in eachindex(range_fc′)
            for idx_as in eachindex(range_as)
                for idx_ec in eachindex(range_dps)
                    for idx_fpe in eachindex(range_fpe)
                        T = range_type[idx_type]
                        fc′ = range_fc′[idx_fc′]
                        fR1 = range_fR1[idx_fc′]
                        fR3 = range_fR3[idx_fc′]

                        as = range_as[idx_as]
                        dps = range_dps[idx_ec]
                        fpe = range_fpe[idx_fpe]
                        #to do
                        #create a Section (ReinforcedConcreteSection or PixelFrameSection <: Section)
                        #by type
                        if T == 3.0
                            compoundsection = make_Y_layup_section(L, t, Lc)
                        elseif T == 2.0
                            compoundsection = make_X2_layup_section(L, t, Lc)
                            #also have to do x4, but will see.
                        elseif T == 4.0
                            compoundsection = make_X4_layup_section(L, t, Lc)
                        else 
                            println("Invalid Type")
                        end
                        # pixelframesection = PixelFrameSection(compoundsection, fc′...) 
                        # pu, mu, vu = get_capacities(pixelframesection)
                        pu, mu, vu, embodied = get_capacities(compoundsection, fc′, fR1, fR3, as, dps, fpe, L)
                        idx_all = [idx_fc′, idx_as, idx_ec, idx_fpe, idx_type]

                        idx = mapping(n, idx_all)
                        results[idx, :] = [fc′, fR1, fR3, as, dps, fpe, pu, mu, vu, embodied, L, t, Lc, T]
                    end
                end
            end
        end
    end

    df = DataFrame(results, [:fc′, :fR1, :fR3, :as, :dps, :fpe, :Pu, :Mu, :Vu, :carbon, :L, :t, :Lc, :T])
    df[!, :ID] = 1:ntotal
    println("Got Catalog with $ntotal outputs")
    return df # results # DataFrame(results)
end

#test
# results_test = get_catalog()
# results = get_catalog(100,10,10,test=false)
# # 11.147s , 575.72 MiB allocation
# date = Dates.today()
# time = Dates.now()

# CSV.write(joinpath(@__DIR__,"Outputs\\output_$date.csv"), results)

results = get_catalog(false)

CSV.write(joinpath(@__DIR__, "Catalogs/FEB1_1_catalog_static.csv"), results)


# calcap(28., 99.0, 0.5, 1600.0)