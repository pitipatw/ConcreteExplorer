include("../Geometry/pixelgeo.jl")
include("../Functions/embodiedCarbon.jl")


"""
Beta 1 value for reinforced concrete compression section calculation
"""
function β1(fc′::Float64)
    out = clamp(0.85 - 0.05 * (fc′ - 28.0) / 7, 0.65, 0.85)
    return out
end

"""
Capacity safety factor ϕ
"""
function ϕ(ϵ::Float64)
    ϕ = clamp(0.65 + 0.25 * (ϵ - 0.002) / 0.003, 0.65, 0.90) # a factor between 0.65 to 0.9
    return ϕ
end

"""
fpy is based on ASTM A421
"""
function get_fps(fpe::Float64, fc′::Float64, ρₚ::Float64; fpy::Float64=1300.0)
    fps = min(fpe + 70.0 + fc′ / (100.0 * ρₚ), fpe + 420.0, fpy)
    return fps
end


""""
get axial capacity of a section
output: P[kN]
"""
function get_Pu(compoundsection::CompoundSection, fc′::Float64, as::Float64, fpe::Float64)

     #concrete area
     ac = compoundsection.area

     # [Alternative] 
     # load section properties file
     # filename = "pixel_$L_$t_$Lc.csv"
     # section = CSV.read(filename, header=true)
 
     #Pure Compression Capacity Calculation
     
     
     ccn = 0.85 * fc′ * ac
     #*need a justification on 0.003 Ep
     # pn = (ccn - (fpe - 0.003 * Ep) * as) / 1000 # convert to [kN]
     pn = (ccn - fpe* as) / 1000 # convert to [kN]
     pu = 0.65 * 0.8 * pn #[kN]
     
     return pu
end

"""
get moment capacity of a post tensioned section
Mu [kNm]
"""
function get_Mu(compoundsection::CompoundSection, fc′::Float64, as::Float64, fpe::Float64, ec::Float64, L::Float64;
    Ep = 200_000,)

    #Pure Moment Capacity
    #concrete area
    ac = compoundsection.area
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

    #rebar position measure from 0.0 (centroid) down, relative value
    rebarpos = ec*(-L)
    #depth is from the top most of the section
    c_depth = depth_from_area(compoundsection,acomp,show_stats = false )
    ymax = compoundsection.ymax #global coordinate

    c_depth_global = ymax - c_depth #global coordinate

    new_sections = Vector{SolidSection}()
    for sub_s in compoundsection.solids
        sub_s_ymax = sub_s.ymax #global coordinate
        sub_s_ymin = sub_s.ymin 
        c_depth_local = sub_s_ymax - c_depth_global
        if c_depth_local > 0
            c_depth_local = clamp(sub_s_ymax - c_depth_global, 0, sub_s_ymax - sub_s_ymin)
            push!(new_sections, sutherland_hodgman(sub_s, c_depth_local, return_section = true))
        end
    end

    cgcomp = CompoundSection(new_sections).centroid

    # depth, cgcomp= getprop(acomp, L, t, Lc)
    
    # clipped_section = sutherland_hodge(section::PolygonalSection, y::Float64; return_section = true)
    # mn_steel = as * fps * arm / 1e6 #[kNm]

    #Recheck with concrete.
    #check compression strain, make sure it's not more than 0.003
    # c = depth
    c = c_depth
    ϵs = fps / Ep
    d = ymax-rebarpos
    ϵc = c * ϵs / ( d - c)

    if ϵc > 0.003
        # Compression strain is more than 0.003
        # recalc based on the compression strain = 0.003
        # ϵs now will be lower than 0.005

        #first find depth based on the 0.003 strain at the top

        ϵc_new = 0.003
        c = L/2 #first guess
        tol = 0.001

        while tol > 0.001
            ϵs_new = ϵc_new*(d - c) / c
            fps_new = ϵs_new * Ep
            acomp = as * fps_new / (0.85 * fc′)
            c_depth = depth_from_area(compoundsection,acomp,show_stats = false )
            tol = abs(c_depth - c)/c
        end
    end



    c_depth_global = ymax - c_depth
    new_sections = Vector{SolidSection}()
    for sub_s in compoundsection.solids
        sub_s_ymax = sub_s.ymax
        sub_s_ymin = sub_s.ymin 

        c_depth_local = sub_s_ymax - c_depth_global
        if c_depth_local > 0
            c_depth_local = clamp(sub_s_ymax - c_depth_global, 0, sub_s_ymax - sub_s_ymin)

                push!(new_sections, sutherland_hodgman(sub_s, c_depth_local, return_section = true))
        end
    end
    
    cgcomp = CompoundSection(new_sections).centroid

    arm = cgcomp[2] - rebarpos
    #moment arm of the section is the distance between the centroid of the compression area and the steel.

    mn = 0.85 * fc′ * ac * arm / 1e6 #[kNm]
    mu = Φ(ϵs) * mn #[kNm]

    return mu
end

"""
find shear capacity based on fiber reinforced
from fib model code.
"""
function get_Vu(compoundsection::CompoundSection, fc′::Float64, as::Float64, fpe::Float64, ec::Float64, L::Float64;
    shear_ratio = 0.30,
    # fR1 = 2.0,
    # fR3 = 2.0 * 0.850)
    fR1 = 3.0,
    fR3 = 3.0) #fR1 and fR3 from PP -> https://eprints.whiterose.ac.uk/179192/3/sustainability-13-11479.pdf

    #Shear calculation.
    ac = compoundsection.area
    d = L
    ashear = ac * shear_ratio
    fctk = 0.17*sqrt(fc′)
    ρs = as / ashear
    k = clamp(sqrt(200.0 / d), 0, 2.0)
    fFts = 0.45 * fR1
    wu = 1.5
    CMOD3 = 1.5
    ptforce = get_Pu(compoundsection, fc′, as, fpe)
    ned = ptforce# can be different
    σcp1 = ned / ac
    σcp2 = 0.2 * fc′
    σcp = clamp(σcp1, 0.0, σcp2)
    fFtu = get_fFtu(fFts, wu, CMOD3, fR1, fR3)
    vn = ashear * get_v(ρs, fc′, fctk, fFtu, 1.0, σcp1, k)# already in kN
    vu = 0.75 * vn
 return vu

end

"""
Bases on Fib Model
"""
function get_fFtu(fFts::Float64, wu::Float64, CMOD3::Float64, fR1::Float64, fR3::Float64)
    fFtu = fFts - wu / CMOD3 * (fFts - 0.5 * fR3 + 0.2 * fR1)
    return fFtu
end

"""
from fib model code
fFtuk = fFtu
fck = fc′
"""
function get_v(ρs::Float64, fck::Float64, fctk::Float64, fFtuk::Float64, γc::Float64, scp::Float64, k::Float64)
    out = 0.18 / γc * k * (100.0 * ρs * (1.0 + 7.5 * fFtuk / fctk) * fck)^(1 / 3) + 0.15 * scp
    return out / 1000.0 #kN
end



"""
Calculate capacities of the given section
Inputs : section information
Outputs:
Pu [kN]
Mu [kNm]
Shear [kN]
"""
function get_capacities(fc′::Float64, as::Float64, ec::Float64, fpe::Float64,
    L::Float64,
    t::Float64,
    Lc::Float64;
    echo = false,
    # L = 102.5, t = 17.5, Lc = 15.,
    # L = 202.5, t = 17.5, Lc = 15.,
    T = "Beam",
    Ep = 200_000,)


    #Calculation starts here.
    
    #Load the right sections (Using AsapSections here)
    if T == "Beam"
        compoundsection = make_Y_layup_section(L, t, Lc)
    elseif T == "Column"
        compoundsection = make_X2_layup_section(L, t, Lc)
        #also have to do x4, but will see.
        # section = make_X4_layup_section(L, t, Lc)
    else
        println("Invalid type")
    end

    # compoundsection = CompoundSection(sections)

    pu = get_Pu(compoundsection, fc′, as, fpe)
    mu = get_Mu(compoundsection, fc′, as, fpe, ec, L)
    vu = get_Vu(compoundsection, fc′, as, fpe, ec, L)

    #Embodied Carbon Calculation
    cfc = fc2e(fc′) #kgCO2e/m3

    # 0.854 kgCo2e per kgsteel
    # 7850 kg/m3
    cst = 0.854*7850 #kgCO2e/m3
    
    ac = compoundsection.area
    embodied = ( ac*cfc + as*cst )/ 1e6 
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


