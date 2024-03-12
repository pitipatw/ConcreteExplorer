include("../Functions/Geometry/pixelgeo.jl")
include("../Functions/embodiedcarbon.jl")

using Printf
"""
Beta 1 value for reinforced concrete compression section calculation
"""
function get_β1(fc′::Float64)
    out = clamp(0.85 - 0.05 * (fc′ - 28.0) / 7, 0.65, 0.85)
    return out
end

"""
Capacity safety factor ϕ
"""
function get_ϕ(ϵ::Float64)
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
function get_Pu(compoundsection::CompoundSection, fc′::Float64, as::Float64, fpe::Float64;
    Ep::Int64 = 200_000)

     #concrete area
     ac = compoundsection.area

     # [Alternative] 
     # load section properties file
     # filename = "pixel_$L_$t_$Lc.csv"
     # section = CSV.read(filename, header=true)
 
     #Pure Compression Capacity Calculation
     
     ccn = 0.85 * fc′ * ac
     #*need a justification on 0.003 Eps
     # pn = (ccn - (fpe - 0.003 * Eps) * as) / 1000 # convert to [kN]
    #  pn = (ccn - (fpe-0.003*Ep)* as)
     #Conservative**, pn will only rely on concrete, not the tendons.
     pn = ccn
     pu = 0.65 * 0.8 * pn /1000 #[kN]
     
     return pu
end

"""
get moment capacity of a post tensioned section
Mu [kNm]
"""
function get_Mu(compoundsection::CompoundSection, fc′::Float64, as::Float64, fpe::Float64, dps::Float64, L::Float64;
            Eps = 200_000,)
    #function get_Mu(pixelFramesection::PixelFrameSection)
    #@unpack compoundsection, fc′, as, fpe, dps, L = pixelframesection
    #as = sum(as)
    #fpe = fpe[1]
    echo = false
    if dps == 0.0
        #work on the centroid.
        #cut the section at centroid and find the maximum moment available from 
        # fc′  = M / Sx. 
        c_global = compoundsection.centroid[2]
        #clip the seciton at depth Y (global)
        # new_sections = Vector{SolidSection}()
        # for sub_s in compoundsection.solids
        #     sub_s_ymax = sub_s.ymax #global coordinate
        #     sub_s_ymin = sub_s.ymin #global coordinate 
        #     if sub_s_ymax - c_global > 0 #this means, c_global is below the top-most point of the subsection.
        #         c_depth_local = clamp(sub_s_ymax - c_global, 0, sub_s_ymax - sub_s_ymin) #it can only go as far as the depth of the section.
        #         push!(new_sections, sutherland_hodgman(sub_s, c_depth_local, return_section = true))
        #     end
        # end
        mn = (fc′- fpe*as/compoundsection.area)*compoundsection.Sx/1e6
        mu = 0.65*mn #compression control section. 
        return mu
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
    acomp = clamp(as * fps / (0.85 * fc′),0, 0.99*ac)

    #recalculating the fps in case it was too large for this section.
    fps = acomp*0.85*fc′/as 

    #depth is from the top most of the section
    compression_depth = depth_from_area(compoundsection,acomp,show_stats = false ) #local
    β1 = get_β1(fc′)
    c_local = compression_depth/β1 #neutral axis

    #find neutral axis at global coordinates
    ymax = compoundsection.ymax #global coordinate max y

    #check concrete compression strain, make sure it's not more than 0.003
    ϵs = fps / Eps
    #rebar position measure from 0.0 (centroid) down, relative value
    rebarpos = -dps  #Centered aroud 0.0, move downward dps [global]
    d = ymax-rebarpos #distance from top of the section to rebar

    if c_local ≥ d 
        """
        An alternative way to analysis this section is required. 
        The current workflow is to bypass it for now. 
        """
        return 0.0
    end
    # section_moment = fps*as*(d-c)

    ϵc = c_local * ϵs / (d - c_local)

    #The strain has to be less than 0.003.
    if ϵc > 0.003
        his_ϵc = ϵc 
        c_local_1 = c_local
        d_1 = d 
        # Compression strain is more than 0.003
        # recalc by forcing the compression strain = 0.003
        # ϵs now will be lower than 0.005
        #first find depth based on the 0.003 strain at the top
        tol = 1
        counter = 0 
        his_local_c = []
        his_ϵs = []
        his_acomp = []
        while tol > 5e-3  #precision problem, if tol is too rough, new ϵc wont be exactly 0.003
            counter  +=1 
            if counter > 1000 || c_local ≥ d 
                # if echo
                println("Counter exceeds limit")
                @show ϵc
                @show fc′, as, dps
                @show c_local_1 
                @show d_1
                @show c_local
                @show d 
                println(his_local_c[end-5:end])
                println(his_ϵs[end-5:end])
                println(his_acomp[end-5:end])
                # end
                return 0.0
                # break
            end
            #numerical value error, let it slightly less than 0.003
            ϵs_new = 0.0029*(d - c_local) / c_local
            fps_new = ϵs_new * Eps
            #limit the compression area between 0 and 0.99 of the actual section.
            acomp = clamp(as * fps_new / (0.85 * fc′),1.0, 0.99*compoundsection.area)
            # @show compoundsection.area
            if acomp <0 
                # if echo
                @show acomp
                @show fps_new
                @show c_local
                @show d
                @show L
                # end
            end
            c_new = depth_from_area(compoundsection,acomp,rel_tol = 1e-3,show_stats = false )
            tol = abs(c_new - c_local)/c_local
            c_local = (c_new+c_local)/2
            ϵs = ϵs_new
            push!(his_local_c,c_local)
            push!(his_ϵs, ϵs)
            push!(his_acomp, acomp)
        end
    end
    #making sure that ϵc is <= 0.003
    ϵc = c_local * ϵs/ (d - c_local)
    @assert ϵc <= 0.003 "$ϵc"

    c_global = ymax - c_local #global coordinate neutral axis

    new_sections = Vector{SolidSection}()
    for sub_s in compoundsection.solids
        sub_s_ymax = sub_s.ymax #global coordinate
        sub_s_ymin = sub_s.ymin #global coordinate 
        if sub_s_ymax - c_global > 0 #this means, c_global is below the top-most point of the subsection.
            c_depth_local = clamp(sub_s_ymax - c_global, 0, sub_s_ymax - sub_s_ymin) #it can only go as far as the depth of the section.
            push!(new_sections, sutherland_hodgman(sub_s, c_depth_local, return_section = true))
        end
    end

    cgcomp = CompoundSection(new_sections).centroid

    arm = cgcomp[2] - rebarpos
    #moment arm of the section is the distance between the centroid of the compression area and the steel.

    mn = 0.85 * fc′ * ac * arm / 1e6 #[kNm]
    mu = get_ϕ(ϵs) * mn #[kNm]

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
    # fFts = 0.45 * fR1
    # wu = 1.5
    # CMOD3 = 1.5
    # fFtuk = fFts - wu / CMOD3 * (fFts - 0.5 * fR3 + 0.2 * fR1)
    fFtuk = fR3/3
    ptforce = fpe*as #get_Pu(compoundsection, fc′, as, fpe)
    ned = 1000*ptforce# N
    σcp1 = ned / ac
    σcp2 = 0.2 * fc′
    σcp = clamp(σcp1, 0.0, σcp2)
 

    #fib notation 
    fck = fc′

    #constant 
    γc = 1.0

    #Fib 2010 gives V design which is Vu. Therefore, no reduction factor of 0.75 is needed.
    vu = ashear * (0.18 / γc * k * (100.0 * ρs * (1.0 + 7.5 * fFtuk / fctk) * fck)^(1 / 3) + 0.15 * σcp)
    
    # vu = 0.75 * vn/1000 # kN
 return vu

end


"""
Calculate capacities of the given section
Inputs : section information
Outputs:
Pu [kN]
Mu [kNm]
Shear [kN]
"""
function get_capacities(compoundsection, fc′::Float64, fR1::Float64, fR3::Float64, as::Float64, dps::Float64, fpe::Float64,
    L::Float64, dosage::Real;
    echo = false)

    pu = get_Pu(compoundsection, fc′, as, fpe)
    mu = get_Mu(compoundsection, fc′, as, fpe, dps, L)
    vu = get_Vu(compoundsection, fc′, fR1, fR3, as, fpe, L,)

    #Embodied Carbon Calculation
    cfc = fc2e(fc′, dosage) #kgCO2e/m3

    # 0.854 kgCo2e per kgsteel
    # 7850 kg/m3
    cst = rebar2e() #0.854*7850 #kgCO2e/m3
    
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


# function get_capacities should give out ϵc. 
