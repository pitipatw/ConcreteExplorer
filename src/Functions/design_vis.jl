using Makie, GLMakie, CairoMakie
include("Geometry/pixelgeo.jl")


function plot_element(element_number::Int64, elements_designs::Dict;
    L::Float64=250.0)

    sections = elements_to_sections[element_number]
    @show ns = length(sections)
    # element_number = string(element_number)

    L = elements_designs[element_number][1][:L]
    t = elements_designs[element_number][1][:t]
    Lc = elements_designs[element_number][1][:Lc]

    @show set_fc′   = [elements_designs[element_number][i][:fc′] for i in 1:ns]
    set_as    = [elements_designs[element_number][i][:as] for i in 1:ns]
    set_fps   = [elements_designs[element_number][i][:fps] for i in 1:ns]
    set_dosage = [elements_designs[element_number][i][:dosage] for i in 1:ns]
    # set_P     = [elements_designs[element_number][i][:load] for i in 1:ns]
    set_axial_force = [elements_designs[element_number][i][:axial_force] for i in 1:ns]
    tendon_profile  = [elements_designs[element_number][i][:dps] for i in 1:ns]
    axial_capacity  = abs.([elements_designs[element_number][i][:Pu] for i in 1:ns])
    moment_capacity = [elements_designs[element_number][i][:Mu] for i in 1:ns]
    shear_capacity  = [elements_designs[element_number][i][:Vu] for i in 1:ns]
    type = [elements_designs[element_number][i][:T] for i in 1:ns]
    if type[1] == 2.0
        type_title = "X2"
    elseif type[1] == 3.0 
        type_title = "Y"
    elseif type[1] == 4.0 
        type_title = "X4"
    else
        print("Invalid type")
    end

    axial_demand  = abs.([demands[i, :pu] for i in sections])
    moment_demand = [demands[i, :mu] for i in sections]
    shear_demand  = [demands[i, :vu] for i in sections]

    #plot center around x = 0 
    n = length(tendon_profile)
    length(axial_demand)
    res = mod(n, 2) * 250
    xmax = div(n, 2) * 500
    x_1 = -xmax + 250 - res
    x_n = xmax - 250 + res
    x_range = x_1:500:x_n

    #create a bands (polygon of possible tendon profile)
    tendon_points = Matrix{Int64}(undef, 2, 2 * n)
    for i in 1:n
        tendon_points[:, i] = [x_range[i], -elements_designs[element_number][i][:min_dps]]
    end
    for i in 1:n
        tendon_points[:, 2*n-i+1] = [x_range[i], -elements_designs[element_number][i][:max_dps]]
    end

    f1 = Figure(size=(1200, 600))
    g = f1[1, 1] = GridLayout()
    axs_design = Axis(g[1, 1], title="Element $type_title $element_number", titlesize=20,
        #aspect=DataAspect(),
        limits=(-1.2*xmax, 1.2*xmax+2*res, -2 * L, 1.2*L+25), #used to be 1.2 L
        yticks=div(L, 100)*100:-100:-div(L, 100)*105,
        # yminorticks = IntervalsBetween(2),
        # yminorgridvisible = true,
        ylabel="y"
    )

    axs_axial = Axis(g[2, 1],# aspect=10,
        limits=(-1.2*xmax, 1.2*xmax+2*res, nothing, nothing), ylabel="Axial [kN]",
    )

    axs_moment = Axis(g[3, 1], #aspect=10,
        limits=(-1.2*xmax, 1.2*xmax+2*res, nothing, nothing), ylabel="Moment [kNm]",
    )
    axs_shear = Axis(g[4, 1], #aspect=10,
        limits=(-1.2*xmax, 1.2*xmax+2*res, nothing, nothing), ylabel="Shear [kN]",
    )

    linkxaxes!(axs_design, axs_axial, axs_moment, axs_shear)

    axs_axial_ratio = Axis(g[2, 1], yaxisposition = :right,# aspect=10,
        limits=(-1.2*xmax, 1.2*xmax+2*res, nothing, nothing), ylabel="Utilization [ratio]",
    )

    axs_moment_ratio = Axis(g[3, 1], yaxisposition =:right, #aspect=10,
        limits=(-1.2*xmax, 1.2*xmax+2*res, nothing, nothing), ylabel="Utilization [ratio]",
    )
    axs_shear_ratio = Axis(g[4, 1], yaxisposition =:right, #aspect=10,
        limits=(-1.2*xmax, 1.2*xmax+2*res, nothing, nothing), ylabel="Utilization [ratio]",
    )

    hidespines!(axs_axial_ratio)
    hidespines!(axs_moment_ratio)
    hidespines!(axs_shear_ratio)
    hideydecorations!(axs_axial_ratio, ticks = false, label =false, ticklabels =false )
    hideydecorations!(axs_moment_ratio, ticks = false, label =false, ticklabels =false )
    hideydecorations!(axs_shear_ratio, ticks = false , label =false, ticklabels =false)
ypos = 0
    points = []
    if type[1] == 3
        for i in 1:n
            points = [x_range[i]-250, -L, 500, L * 1.5]
            poly!(axs_design, Rect(points...), color=set_fc′[i], 
            colorrange= (0,80), colormap=cgrad(:Greys_3), 
            strokecolor= :black,strokewidth = 2)
            text!(axs_design,points[1], points[2], text = string(set_fc′[i]))
            text!(axs_design,points[1], points[2]+50, text = string(set_dosage[i]))

        end
        section_pts = make_Y_layup_section(L,t,Lc)
        for s in section_pts.solids
            poly!(axs_design, s.points.+[x_n+res+500;0], color = :grey)
        end
        ypos = 150
    elseif type[1] == 2 
        for i in 1:n
            points = [x_range[i]-250, -0.5*L, 500, L]
            poly!(axs_design, Rect(points...), color=set_fc′[i], 
            colorrange= (0,80), colormap=cgrad(:Greys_3), 
            strokecolor= :black,strokewidth = 2)
            text!(axs_design,points[1], points[2], text = string(set_fc′[i]))
            text!(axs_design,points[1], points[2]+50, text = string(set_dosage[i]))

        end
        #plot cross section on the side of the beam plot
        section_pts = make_X2_layup_section(L,t,Lc)
        for s in section_pts.solids
            poly!(axs_design, s.points.+[x_n+res+500;0], color = :grey)
            
        end
        ypos = -250
    elseif type[1] == 4
        for i in 1:n
            points = [x_range[i]-250, -L, 500, 2*L]
            poly!(axs_design, Rect(points...), color=set_fc′[i], 
            colorrange= (0,80), colormap=cgrad(:Greys_3), 
            strokecolor= :black,strokewidth = 2)
            text!(axs_design,points[1], points[2], text = string(set_fc′[i]))
            text!(axs_design,points[1], points[2]+50, text = string(set_dosage[i]))

        end
        #plot cross section on the side of the beam plot
        section_pts = make_X4_layup_section(L,t,Lc)
        for s in section_pts.solids
            poly!(axs_design, s.points.+[x_n+res+500;0], color = :grey)
        end
        ypos = -350
    else
        println("Invalid Type")
    end
        
    Colorbar(f1[1, 2], ticks=[40,60,80], colorrange=(40, 100), colormap=cgrad(:Greys_3))


    tendon = poly!(axs_design, tendon_points, color=:skyblue, alpha=0.5, transparent=true)
    tendon_pts = scatter!(axs_design, tendon_points)
    tendon = lines!(axs_design, x_range, -tendon_profile, linewidth = 5)

    #add text description
    txt_force = round(maximum(set_fps.*set_as)/1000,digits = 3) #kN.
    
    text!(axs_design,-xmax-500, ypos, text = "Tendon load: "*string(txt_force)*" kN")
    # text!(axs_design,-xmax+700, ypos, text = "Apply test load "*string(set_P[1])* " kN")
    text!(axs_design,-xmax+2000,ypos, text = "Axial post tensioned force "*string(round(set_axial_force[1], digits = 3))* " kN")
    

    # hidexdecorations!(axs_design, grid=false)
    # hidexdecorations!(axs_axial, grid=false)
    # hidexdecorations!(axs_moment, grid=false)

    # lines!(axs_axial, x_range, axial_capacity, color=:red)
    # lines!(axs_moment, x_range, moment_capacity, color=:blue)
    # lines!(axs_shear, x_range, shear_capacity, color=:green)

    # lines!(axs_axial, x_range, axial_demand, linestyle=:dash, color=:red)
    # lines!(axs_moment, x_range, moment_demand, linestyle=:dash, color=:blue)
    # lines!(axs_shear, x_range, shear_demand, linestyle=:dash, color=:green)

    stairs!(axs_axial, x_range, axial_capacity, color=:red, step=:center)
    stairs!(axs_moment, x_range, moment_capacity, color=:blue, step=:center)
    stairs!(axs_shear, x_range, shear_capacity, color=:green, step=:center)

    stairs!(axs_axial, x_range, axial_demand, linestyle=:dash, color=:red, step=:center)
    stairs!(axs_moment, x_range, moment_demand, linestyle=:dash, color=:blue, step=:center)
    stairs!(axs_shear, x_range, shear_demand, linestyle=:dash, color=:green, step=:center)

    # scatter!(axs_axial, x_range, axial_capacity, color=:red, step=:center)
    # scatter!(axs_moment, x_range, moment_capacity, color=:blue, step=:center)
    # scatter!(axs_shear, x_range, shear_capacity, color=:green, step=:center)

    # scatter!(axs_axial, x_range, axial_demand, color=:red, step=:center, marker = 'x', markersize = 20)
    # scatter!(axs_moment, x_range, moment_demand, color=:blue, step=:center, marker = 'x', markersize = 20)
    # scatter!(axs_shear, x_range, shear_demand, color=:green, step=:center, marker = 'x', markersize = 20)


    axial_utilization  = [ axial_demand[i] ≤ 0.1 ? 0 : axial_capacity[i]/axial_demand[i] for i in eachindex(axial_demand)]
    moment_utilization = [ moment_demand[i] ≤ 0.1 ? 0 : moment_capacity[i]/moment_demand[i] for i in eachindex(moment_demand)]
    shear_utilization  = [ shear_demand[i] ≤ 0.1 ? 0 : shear_capacity[i]/shear_demand[i] for i in eachindex(shear_demand)]

    # moment_utilization =
    # shear_utilizatishear
    scatter!(axs_axial_ratio, x_range, axial_utilization, color=:red, marker = :utriangle)
    scatter!(axs_moment_ratio, x_range, moment_utilization, color=:blue, marker = :utriangle)
    scatter!(axs_shear_ratio, x_range, shear_utilization, color=:green, marker = :utriangle)



    # for (i, label) in enumerate(["Axial [kN]", "Moment [kNm]", "Shear [kN]"])
    #     Box(g[i, 2], color = :gray90)
    #     Label(g[i,2], label, rotation = pi/2, tellheight = false)
    # end

    rowgap!(g, 10)

    yspace = maximum(tight_yticklabel_spacing!, [axs_axial, axs_shear, axs_moment]) + 15

    axs_axial.yticklabelspace = yspace
    axs_moment.yticklabelspace = yspace
    axs_shear.yticklabelspace = yspace



    return f1
end