using AsapSections
using Makie , GLMakie
include("../Geometry/pixelgeo.jl")
include("postTensionedFunc.jl")

# function testcapacity()
    #inputs
    fc′ = 25.0
    as = 200.0
    ec = 0.5
    fpe = 1200.0
    L = 200.0
    t = 15.0
    Lc = 15.0
    Ep = 200_000

    compoundsection = make_Y_layup_section(L, t, Lc)
    # compoundsection = c1
    f1 = Figure(resolution = (500,500))
    ax1 = Axis(f1[1,1], aspect = DataAspect(), title = "test plot")
    for s in compoundsection.solids
        println(s.points)
        lines!(ax1, s.points)
    end

    ac = compoundsection.area

    #Pure Compression Capacity Calculation

    ccn = 0.85 * fc′ * ac
    #*need a justification on 0.003 Ep
    pn = (ccn - (fpe - 0.003 * Ep) * as) / 1000 # convert to [kN]
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
    areas = Vector{Float64}()
    depths = Vector{Float64}()
    # for i in 100:100:compoundsection.area
    # acomp = 3000.0
    # acomp =i
    if acomp >= ac
        println("Acomp exceeds Ac, using Ac instead")
        acomp = ac
        c_depth =  compoundsection.ymax- compoundsection.ymin
    else
        c_depth = depth_from_area(compoundsection, acomp, show_stats=false)
    end

    #rebar position measure from 0.0 (centroid) down, relative value
    rebarpos = ec * (-L)
    #depth is from the top most of the section

# f2 = Figure(resolution = (200,200))
# ax10 = Axis(f2[1,1])
# scatter!(ax10, depths,areas)
    # depth, cgcomp= getprop(acomp, L, t, Lc)

    # clipped_section = sutherland_hodge(section::PolygonalSection, y::Float64; return_section = true)
    # mn_steel = as * fps * arm / 1e6 #[kNm]

    #Recheck with concrete.
    #check compression strain, make sure it's not more than 0.003
    # c = depth
    ymax = compoundsection.ymax #global coordinate

    c = c_depth
    ϵs = fps / Ep
    d = ymax - rebarpos
    ϵc = c * ϵs / (d - c)

    cached_ϵc = ϵc


    if ϵc > 0.003
        # Compression strain is more than 0.003
        # recalc based on the compression strain = 0.003
        # ϵs now will be lower than 0.005

        #first find depth based on the 0.003 strain at the top
        ϵc_new = 0.003

        c = 0.003*d/(ϵs+0.003)
        tol = 1.0
        while tol > 0.01  #limited by Asap's tol of 0.001
            ϵs_new = 0.003* (d - c) / c
            fps_new = ϵs_new * Ep
            acomp = as * fps_new / (0.85 * fc′)
            if acomp > compoundsection.area
                acomp = compoundsection.area          
            end
            c_new = depth_from_area(compoundsection, acomp, show_stats=true)
            @show tol = abs(c_new - c) / c
            @show c = (c_new+c)/2
            # @show c = 0.003*d/(ϵs_new+0.003)
            # @show  ϵc_new = c * ϵs / (d - c)

        end
    end

    c
    ϵs_new = 0.003* (d - c) / c
    ϵc = c * ϵs_new / (d - c)


    #after we got the correct compression depth, we get cg of that compression area.
    c_depth_global = ymax - c #global coordinate
    # hlines!(ax1, c_depth_global)

    ax2 = Axis(f1[2,1], aspect = DataAspect(), limits = ax1.limits)
    # ax3 = Axis(f1[2,2], aspect = DataAspect(), limits = ax1.limits)
    # ax4 = Axis(f1[2,3],aspect = DataAspect(), limits = ax1.limits)
    # axs = [ax2, a]

    new_sections = Vector{SolidSection}()
    for k in 1:length(compoundsection.solids)
        sub_s = compoundsection.solids[k]
        sub_s_ymax = sub_s.ymax  #global coordinate
        sub_s_ymin = sub_s.ymin
        maxd = sub_s_ymax - sub_s_ymin
        c_depth_local = sub_s_ymax - c_depth_global
        if sub_s_ymax >= c_depth_global
            if c_depth_global < sub_s_ymin
                c_depth_local = sub_s_ymax - sub_s_ymin
            end
            # @show k
            new_sub_sec = sutherland_hodgman(sub_s, c_depth_local, return_section = true)
            # @show length(new_sub_sec.points)
            # if k == 2 
                # println(new_sub_sec.points)
            # println(new_sub_sec.area)
            lines!(ax2, new_sub_sec.points)
            push!(new_sections, sutherland_hodgman(sub_s, c_depth_local, return_section=true))
        end
    end
    f1
    
    clipped_compoundsection = CompoundSection(new_sections)
    
    
    # push!(areas,clipped_compoundsection.area)
    # push!(depths ,c_depth_global )
    cgcomp = clipped_compoundsection.centroid
    # scatter!(ax1, Point2(cgcomp), color = :red)

    arm = cgcomp[2] - rebarpos
    #moment arm of the section is the distance between the centroid of the compression area and the steel.

    mn = 0.85 * fc′ * ac * arm / 1e6 #[kNm]
    mu = Φ(ϵs) * mn #[kNm]


# end

# testcapacity()