using Dates
using DataFrames
using CSV
# using BenchmarkTools
using AsapToolkit
using Printf

include("Geometry/pixelgeo.jl")
include("embodiedCarbon.jl")
include("capacities.jl")

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

"""
    get_catalog(test::Bool)::DataFrame
"""
function get_catalog(case::String; 
    fc′_path::String = "src//Tables//fiber_with_extrapolation.csv")::DataFrame

    println("Getting fc′ info from ", fc′_path)
    pixel_sections = [205.0 35.0 30.0;] #updated by Jenna Jan 2024
    nps = size(pixel_sections)[1]
    println("Getting catalog with $nps section(s)")
    out = DataFrame()
    for i in 1:size(pixel_sections)[1]
        L, t, Lc = pixel_sections[i, :]
        sub_catalog = get_catalog(L, t, Lc, case = case, fc′_path = fc′_path)
        out = vcat(out, sub_catalog)
    end
    return out
end

"""
    get_catalog(L,t,Lc; test = true)::DataFrame
test = true (depreciated) get a dummy catalog.
test = false get a full catalog.
"""
function get_catalog(L, t, Lc;
     case::String = "test", fc′_path::String = "src//Tables//fiber_with_extrapolation.csv")::DataFrame
    if case == "test" #Depreciated
        #test
        println("Running test case")
        range_fc′ = 28.0
        range_as = 140.0
        range_dps = 50.0
        range_fpe = 186.0
    elseif case == "default" 
        println("Running default mode")
        fc_fiber = CSV.read(fc′_path, DataFrame)
        #These come in a set of (fc′,fR1, fR3)
        range_fc′= convert(Array{Float64}, fc_fiber[!, :strength]) #some numbers can be Int.
        range_fR1 = fc_fiber[!, :fR1]
        range_fR3 = fc_fiber[!, :fR3]
        range_dosage = fc_fiber[!, :dosage]
        @assert length(range_fc′) == length(range_fR1) "Error! Number of rows of fc′ ≠ number of rows of fR1 "
        @assert length(range_fc′) == length(range_fR3) "Error! Number of rows of fc′ ≠ number of rows of fR3 "
        @assert length(range_fc′) == length(range_dosage) "Error! Number of rows of fc′ ≠ number of rows of fiber dosage "

        range_as = [(99.0 * 2),(140.0 * 2)] # x2 are for 2 ropes on 2 sides 12.7 and 15.2 mm dia wires.

        range_dps_2 = [0.0]
        range_dps_3 = vcat(0.0:20:350.0) 
        range_dps_4 = [0.0]

        range_fpe = (0.00:0.050:0.5) * 1860.0 #MPa
        range_type = [2.0, 3.0, 4.0] #PixelFrame configuration -> Y = 3 ,X2 = 2, X4 = 4.
    
    else
        println("Error Invalid test case")
        return nothing
    end

    # n = length.([range_fc′, range_as, range_dps, range_fpe, range_type])
    # @show ntotal = prod(n)
    #Pre allocating results
    @show n_result2 = length.([range_dps_2,range_fc′, range_as,  range_fpe])
    @show n_result3 = length.([range_dps_3,range_fc′, range_as,  range_fpe])
    @show n_result4 = length.([range_dps_4,range_fc′, range_as,  range_fpe])
    @show n_total = sum(prod.([n_result2, n_result3, n_result4]))
    results = Matrix{Float64}(undef, n_total, 10 + 2 + 3) #L t Lc
    #we will loop through these three parameters and get the results.
    # with constant cross section properties.
    for idx_type in eachindex(range_type)
        T = range_type[idx_type]
        if T == 3 #Y layup 
            range_dps = range_dps_3
            compoundsection = make_Y_layup_section(L, t, Lc)
        elseif T == 2.0 # X2 layup
            range_dps = range_dps_2
            compoundsection = make_X2_layup_section(L, t, Lc)
        elseif T == 4.0 # X4 layup
            range_dps = range_dps_4
            compoundsection = make_X4_layup_section(L, t, Lc)
        else
            println("Warning, Type not supported")
        end
        for idx_dps in eachindex(range_dps)
            for idx_fc′ in eachindex(range_fc′)
                for idx_as in eachindex(range_as)
                    for idx_fpe in eachindex(range_fpe)

                        fc′ = range_fc′[idx_fc′]
                        fR1 = range_fR1[idx_fc′]
                        fR3 = range_fR3[idx_fc′]
                        dosage = range_dosage[idx_fc′]

                        as = range_as[idx_as]
                        dps = range_dps[idx_dps]
                        fpe = range_fpe[idx_fpe]
                        #to do
                        #create a Section (ReinforcedConcreteSection or PixelFrameSection <: Section)
                        # pixelframesection = PixelFrameSection(compoundsection, fc′...) 
                        # pu, mu, vu = get_capacities(pixelframesection)
                        pu, mu, vu, embodied = get_capacities(compoundsection, fc′, fR1, fR3, as, dps, fpe, L, dosage)


                        idx_all = [idx_dps, idx_fc′, idx_as, idx_fpe]

                        if T == 2 
                            sub_n = n_result2
                            shift = 0
                        elseif T == 3 
                            sub_n = n_result3
                            shift = prod(n_result2)
                        elseif T == 4
                            sub_n = n_result4
                            shift = prod(n_result2) + prod(n_result3)
                        end
                        idx = mapping(sub_n, idx_all)
                        # @show idx, T
                        results[idx + shift, :] = [fc′,dosage, fR1, fR3, as, dps, fpe, pu, mu, vu, embodied, L, t, Lc, T]
                    end
                end
            end
        end
    end

    df = DataFrame(results, [:fc′, :dosage,:fR1, :fR3, :as, :dps, :fpe, :Pu, :Mu, :Vu, :carbon, :L, :t, :Lc, :T])
    df[!, :ID] = 1:n_total
    println("Got Catalog with $n_total outputs")
    return df # results # DataFrame(results)
end




