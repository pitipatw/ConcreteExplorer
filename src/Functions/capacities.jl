include("../Functions/Geometry/pixelgeo.jl")
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
     pn = (ccn - fpe* as)
     pu = 0.65 * 0.8 * pn /1000 #[kN]
     
     return pu
end

"""
get moment capacity of a post tensioned section
Mu [kNm]
"""
function get_Mu(compoundsection::CompoundSection, fc′::Float64, as::Float64, fpe::Float64, dps::Float64, L::Float64;
    Ep = 200_000,)
#function get_Mu(pixelFramesection::PixelFrameSection)
#@unpack compoundsection, fc′, as, fpe, dps, L = pixelframesection
#as = sum(as)
#fpe = fpe[1]
    if dps == 0.0
        return 0.0
    end 
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
        #re-calculate fps 
        fps = acomp*0.85*fc′/as 
        @assert acomp = as * fps / (0.85 * fc′)
    end


    #depth is from the top most of the section
    c = depth_from_area(compoundsection,acomp,show_stats = false ) #local

    #find depth at global coordinates
    ymax = compoundsection.ymax #global coordinate
    c_depth_global = ymax - c #global coordinate

    new_sections = Vector{SolidSection}()
    for sub_s in compoundsection.solids
        sub_s_ymax = sub_s.ymax #global coordinate
        sub_s_ymin = sub_s.ymin #global coordinate
        c_depth_local = sub_s_ymax - c_depth_global
        if c_depth_local > 0 #this means, c_depth_global is below the top part of the subsection.
            c_depth_local = clamp(sub_s_ymax - c_depth_global, 0, sub_s_ymax - sub_s_ymin) #it can only go as far as the depth of the section.
            push!(new_sections, sutherland_hodgman(sub_s, c_depth_local, return_section = true))
        end
    end

    #Recheck with concrete.
    #check compression strain, make sure it's not more than 0.003
    ϵs = fps / Ep
    #rebar position measure from 0.0 (centroid) down, relative value
    rebarpos = 0.0 - dps  #Centroud at 0.0, move downward dps
    d = ymax-rebarpos
    @assert d > 0
    @assert d > c
    ϵc = c * ϵs / (d - c)

    if ϵc > 0.003
        # Compression strain is more than 0.003
        # recalc by forcing the compression strain = 0.003
        # ϵs now will be lower than 0.005

        #first find depth based on the 0.003 strain at the top

        tol = 1
        counter = 0 
        while tol > 1e-3  #precision problem, if tol is too rough, new ϵc wont be exactly 0.003
            counter  +=1 
            if counter > 1000
                println("Counter exceeds limit")
                @show fc′, as, dps
                return 0.0
                break
            end
            #numerical value error, let it slightly less than 0.003
            ϵs_new = 0.0029*(d - c) / c
            fps_new = ϵs_new * Ep
            #limit the compression area between 0 and 0.99 of the actual section.
            acomp = clamp(as * fps_new / (0.85 * fc′),1.0, 0.99*compoundsection.area)
            # @show compoundsection.area
            if acomp <0 
                @show acomp
                @show fps_new

                @show c
                @show d
                @show L

            end
            c_new = depth_from_area(compoundsection,acomp,rel_tol = 1e-3,show_stats = false )
            tol = abs(c_new - c)/c
            c = c_new
            ϵs = ϵs_new
        end
        #making sure that ϵc is <= 0.003
        ϵc = c * ϵs/ (d - c)
        # @show ϵc
        @assert ϵc <= 0.003

        c_depth_global = ymax - c
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

        #check es again to see if es can 
        # @show ϵs
        @assert ϵs >= 0.002

    end
    
    cgcomp = CompoundSection(new_sections).centroid

    # @show rebarpos
    # @show cgcomp
    arm = cgcomp[2] - rebarpos
    #moment arm of the section is the distance between the centroid of the compression area and the steel.

    mn = 0.85 * fc′ * ac * arm / 1e6 #[kNm]
    mu = ϕ(ϵs) * mn #[kNm]

    return mu
end

"""
find shear capacity based on fiber reinforced
from fib model code.
"""
function get_Vu(compoundsection::CompoundSection, fc′::Float64, fR1::Float64, fR3::Float64, as::Float64, fpe::Float64, L::Float64;
    shear_ratio = 1.0)
    # fR1 = 2.0,
    # fR3 = 2.0 * 0.850)
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
    ned = 1000*ptforce# can be different
    σcp1 = ned / ac
    σcp2 = 0.2 * fc′
    σcp = clamp(σcp1, 0.0, σcp2)
    fFtuk = fFts - wu / CMOD3 * (fFts - 0.5 * fR3 + 0.2 * fR1)

    #fib notation 
    fck = fc′

    #constant 
    γc = 1.0

    vn = ashear * (0.18 / γc * k * (100.0 * ρs * (1.0 + 7.5 * fFtuk / fctk) * fck)^(1 / 3) + 0.15 * σcp)
    
    vu = 0.75 * vn/1000 # kN

 return vu

end
#get Vn , then Vu = 0.75Vn

# function get_Vu(compoundsection::CompoundSection, fc′::Float64, fR1::Float64, fR3::Float64, as::Float64, fpe::Float64, dps::Float64, L::Float64)

#     a = 1 
#     return
# end



"""
Calculate capacities of the given section
Inputs : section information
Outputs:
Pu [kN]
Mu [kNm]
Shear [kN]
"""
function get_capacities(compoundsection, fc′::Float64, fR1::Float64, fR3::Float64, as::Float64, dps::Float64, fpe::Float64,
    L::Float64;
    echo = false)

    #Calculation starts here.
    
    # #Load the right sections (Using AsapSections here)
    # if T == "Beam"
    #     compoundsection = make_Y_layup_section(L, t, Lc)
    # elseif T == "Column"
    #     compoundsection = make_X2_layup_section(L, t, Lc)
    #     #also have to do x4, but will see.
    #     # section = make_X4_layup_section(L, t, Lc)
    # else
    #     println("Invalid type")
    # end

    # compoundsection = CompoundSection(sections)

    pu = get_Pu(compoundsection, fc′, as, fpe)
    mu = get_Mu(compoundsection, fc′, as, fpe, dps, L)
    vu = get_Vu(compoundsection, fc′, fR1, fR3, as, fpe, L,)

    #Embodied Carbon Calculation
    cfc = fc2e(fc′) #kgCO2e/m3

    # 0.854 kgCo2e per kgsteel
    # 7850 kg/m3
    cst = 0.854*7850 #kgCO2e/m3
    
    ac = compoundsection.area
    embodied = ( ac*cfc + as*cst )/ 1e6 # mm2 -> m2
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


