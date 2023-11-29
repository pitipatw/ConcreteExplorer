# using Dates
# using BenchmarkTools
using AsapSections
using Printf
include("Geometry/pixelgeo.jl")


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
Map n dimentional vector into an index.
"""
function mapping(n::Vector{Int64}, idx::Vector{Int64})
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
Get embodied carbon coefficient of concrete based on fc′
input : fc′ [MPa]
output: ecc of fc′ [kgCO2e/m3]
"""
function fc2e(fc′::Real)
    out = -0.0626944435544512 * fc′^2 + 10.0086510099949 * fc′ + 84.14807
    return out
end

"""
Calculate capacities of the given section
Inputs : section information
Outputs:
Pu [kN]
Mu [kNm]
Shear [kN]
"""
function get_capacities(fc′, as, ec, fpe, L, t, Lc;
    echo=false,
    # L = 102.5, t = 17.5, Lc = 15.,
    # L = 202.5, t = 17.5, Lc = 15.,
    T="Beam",
    Ep=200_000,
    shear_ratio=0.30,
    fR1=2.0,
    fR3=2.0 * 0.850)

    #Calculation starts here.

    if T == "Beam"
        section = make_Y_layup_section(L, t, Lc)

        ac = section.area
        # y, A = depth_map(compoundsection, 250)

    elseif T == "Column"
        section = make_X2_layup_section(L, t, Lc)

        #also have to do x4, but will see how can we do this, run 2 times and embeded?
        # section = make_X4_layup_section(L, t, Lc)
        # compoundsection = CompoundSection(section)
        ac = section.area
        # y, A = depth_map(compoundsection, 250)
    else
        println("Invalid type")
    end

    #load section properties file
    # filename = "pixel_$L_$t_$Lc.csv"
    # section = CSV.read(filename, header=true)


    #Pure Compression Capacity
    ccn = 0.85 * fc′ * ac
    #need a justification on 0.003 Ep
    pn = (ccn - (fpe - 0.003 * Ep) * as) / 1000 #[kN]
    pu = 0.65 * 0.8 * pn #[kN]
    ptforce = pu #[kN]


    #Pure Moment Capacity

    #From ACI318M-19 Table: 20.3.2.4.1
    ρ = as / ac #reinforcement ratio (Asteel/Aconcrete)
    fps1 = fpe + 70 + fc′ / (100 * ρ) #
    fps2 = fpe + 420
    fps3 = 1300.0 #Yield str of steel from ASTM A421
    fps = minimum([fps1, fps2, fps3])

    #concrete compression area balanced with steel tension force.
    acomp = as * fps / (0.85 * fc′)
    if acomp > ac
        println("Acomp exceeds Ac, using Ac instead")
        acomp = ac
    end
    #rebar position measure from 0.0 down
    rebarpos = ec * (-L)
    @show d = depth_from_area(section, acomp)

    # depth, cgcomp= getprop(acomp, L, t, Lc)

    # clipped_section = sutherland_hodge(section::PolygonalSection, y::Float64; return_section = true)
    # mn_steel = as * fps * arm / 1e6 #[kNm]

    #Recheck with concrete.
    #check compression strain, make sure it's not more than 0.003
    # c = depth
    c = d
    ϵs = fps / Ep
    ϵc = c * ϵs / (rebarpos - c)

    if ϵc > 0.003
        # println("Compression strain is $ϵc which is more than 0.003")
        #recalc based on the compression strain
        #calculate the new compression area based on the strain
        #first find depth based on the 0.003 strain at the top

        ϵc_new = 0.003
        depth = L #first guess
        tol = 0.001

        while tol > 0.001
            ϵs_new = ϵc_new * (rebarpos - depth) / depth
            fps_new = ϵs_new * Ep

            acomp = as * fps_new / (0.85 * fc′)

            # depth_new, cgcomp= getprop(acomp, L, t, Lc)
            tol = abs(depth_new - depth) / depth
            depth = depth_new
        end
    end
    cgcomp = 100 #dummy
    arm = cgcomp - rebarpos
    #moment arm of the section is the distance between the centroid of the compression area and the steel.

    mn = 0.85 * fc′ * ac * arm / 1e6 #[kNm]
    mu = Φ(ϵs) * mn #[kNm]

    # println("#"^50)

    #Shear Calculation
    d = L
    ashear = ac * shear_ratio
    fctk = 0.17 * sqrt(fc′)
    ρs = as / ashear
    k = clamp(sqrt(200.0 / d), 0, 2.0)
    fFts = 0.45 * fR1
    wu = 1.5
    CMOD3 = 1.5
    ned = ptforce# can be different
    σcp1 = ned / ac
    σcp2 = 0.2 * fc′
    σcp = clamp(σcp1, 0.0, σcp2)
    fFtu = get_fFtu(fFts, wu, CMOD3, fR1, fR3)
    vn = ashear * get_v(ρs, fc′, fctk, fFtu, 1.0, σcp1, k)# already in kN
    vu = 0.75 * vn


    #Embodied Carbon Calculation
    cfc = fc2e(fc′) #kgCO2e/m3
    # 0.854 #kgper kg
    # 7850 kg/m3
    cst = 0.854 * 7850 #kgCO2e/m3

    embodied = (ac * cfc + as * cst) / 1e6
    if echo
        @printf "The pure compression capacity is %.3f [kN]\n" pu
        @printf "The pure moment capacity is %.3f [kNm]\n" mu
        @printf "The shear capacity is %.3f [kN]\n" vu
        @printf "The embodied carbon is %.3f [kgCo2e/m3]" embodied
    end

    #write output into CSV
    # dataall = hcat(val,res,checkres)
    # table1 = DataFrame(dataall, :auto)
    # CSV.write("output.csv", table1)

    #parallel plot the result 

    #Scatter plot the result.
    return pu, mu, vu, embodied
end


function get_catalog()
    L = 102.5
    t = 17.5
    Lc = 15.0
    return get_catalog(L, t, Lc)
end

function get_catalog(L, t, Lc; run_test=true)
    if !run_test
        range_fc′ = 28.0:7.0:56.0
        range_as = [99.0, 140.0]
        range_ec = 0.5:0.1:1.2
        range_fpe = (0.00:0.02:0.5) * 1860.0
    else
        #test
        range_fc′ = 28.0
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

                    idx = mapping(n, idx_all)
                    results[idx, :] = [fc′, as, ec, fpe, pu, mu, vu, embodied]
                end
            end
        end
    end
    df = DataFrame(results, [:fc′, :as, :ec, :fpe, :Pu, :Mu, :Vu, :carbon])
    df[!, :ID] = 1:total_s
    return df# results # DataFrame(results)
end

#test
results_test = get_catalog()
results = get_catalog(100, 10, 10, run_test=false)
# 11.147s , 575.72 MiB allocation
date = Dates.today()
time = Dates.now()

# CSV.write(joinpath(@__DIR__, "Outputs\\output_$date.csv"), results)


# calcap(28., 99.0, 0.5, 1600.0)