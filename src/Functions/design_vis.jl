using Makie, GLMakie, CairoMakie
include("Geometry/pixelgeo.jl")

function plot_feasible_sections(all_feasible_sections, save_image_directory)
    # Threads.nthreads()
    if !ispath(save_image_directory * "check_sections//")
        mkdir(save_image_directory * "check_sections//")
    end
    figures = Vector{Figure}(undef, size(demands)[1])
    #select a section to see the available designs
    # Threads.@threads for section_number in 1:size(demands)[1]
    for section_number in 1:size(demands)[1]
        if section_number < 100
            @show section_number
        end
        # section_number = 1

        figure_check_section = Figure(size=(500, 500))
        ax_1 = Axis(figure_check_section[1, 1], xlabel="Moment [kNm]", ylabel="Shear [kN]", title=string(section_number))
        ax_2 = Axis(figure_check_section[2, 1], xlabel="dps", ylabel="T")

        scatter!(ax_1, catalog[all_feasible_sections[section_number], :Mu], catalog[all_feasible_sections[section_number], :Vu], color=catalog[all_feasible_sections[section_number], :fc′], colorrange=extrema(catalog[!, :fc′]))
        scatter!(ax_2, catalog[all_feasible_sections[section_number], :dps], catalog[all_feasible_sections[section_number], :T], color=catalog[all_feasible_sections[section_number], :fc′], colorrange=extrema(catalog[!, :fc′]))

        scatter!(ax_1, demands[section_number, :mu], demands[section_number, :vu], marker='x')
        type_map = Dict("primary" => 3, "secondary" => 3, "columns" => 2)
        scatter!(ax_2, demands[section_number, :ec_max]*1000, getindex.(Ref(type_map), demands[section_number, :type]), marker='x')

        figures[section_number] = figure_check_section
        # figure_check_section
        # @show imagesavepath * "figure_check_section" * string(section_number) * ".png"
        # save(imagesavepath * "figure_check_section" * string(section_number) * ".png", figure_check_section)
    end

    for i in eachindex(figures)
        save(save_image_directory * "check_sections//figure_check_section_" * string(i) * ".png", figures[i])
    end

    for i in eachindex(all_feasible_sections)
        print("Section $i")
        println(length(all_feasible_sections[i]))
    end

    return
end


function plot_element(element_number::Int64, elements_designs::Dict, elements_to_sections;
    L::Float64=250.0)

    sections = elements_to_sections[element_number]
    ns = length(sections)
    # element_number = string(element_number)

    L = elements_designs[element_number][1][:L]
    t = elements_designs[element_number][1][:t]
    Lc = elements_designs[element_number][1][:Lc]

    @show set_fc′ = [elements_designs[element_number][i][:fc′] for i in 1:ns]
    set_as = [elements_designs[element_number][i][:as] for i in 1:ns]
    set_fps = [elements_designs[element_number][i][:fps] for i in 1:ns]
    set_dosage = [elements_designs[element_number][i][:dosage] for i in 1:ns]
    # set_P     = [elements_designs[element_number][i][:load] for i in 1:ns]
    set_axial_force = [elements_designs[element_number][i][:axial_force] for i in 1:ns]
    tendon_profile = [elements_designs[element_number][i][:dps] for i in 1:ns]
    axial_capacity = abs.([elements_designs[element_number][i][:Pu] for i in 1:ns])
    moment_capacity = [elements_designs[element_number][i][:Mu] for i in 1:ns]
    shear_capacity = [elements_designs[element_number][i][:Vu] for i in 1:ns]
    set_embodied_carbon = [elements_designs[element_number][i][:carbon] for i in 1:ns]
    type = [elements_designs[element_number][i][:T] for i in 1:ns]

    set_max_dps = 1000 .*[demands[elements_to_sections[element_number][i], :ec_max] for i in eachindex(elements_to_sections[ element_number])]
    if type[1] == 2.0
        type_title = "X2"
    elseif type[1] == 3.0
        type_title = "Y"
    elseif type[1] == 4.0
        type_title = "X4"
    else
        print("Invalid type")
    end

    axial_demand = abs.([demands[i, :pu] for i in sections])
    moment_demand = [demands[i, :mu] for i in sections]
    shear_demand = [demands[i, :vu] for i in sections]

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
    axs_design = Axis(g[1, 1], title="Element $type_title $element_number (L,t,Lc): ($L,$t,$Lc)", titlesize=20,
        #aspect=DataAspect(),
        # limits=(-1.2*xmax, 1.2*xmax+2*res, -2 * L, 1.2*L+25), #used to be 1.2 L
        limits=(-1.2 * xmax, 1.2 * xmax + 2 * res, -400, 1.2 * L + 100), #used to be 1.2 L
        yticks=div(L, 100)*100:-100:-div(L, 100)*105,
        # yminorticks = IntervalsBetween(2),
        # yminorgridvisible = true,
        ylabel="y"
    )

    axs_axial = Axis(g[2, 1],# aspect=10,
        limits=(-1.2 * xmax, 1.2 * xmax + 2 * res, nothing, nothing), ylabel="Axial [kN]",
    )

    axs_moment = Axis(g[3, 1], #aspect=10,
        limits=(-1.2 * xmax, 1.2 * xmax + 2 * res, nothing, nothing), ylabel="Moment [kNm]",
    )
    axs_shear = Axis(g[4, 1], #aspect=10,
        limits=(-1.2 * xmax, 1.2 * xmax + 2 * res, nothing, nothing), ylabel="Shear [kN]",
    )

    linkxaxes!(axs_design, axs_axial, axs_moment, axs_shear)

    axs_axial_ratio = Axis(g[2, 1], yaxisposition=:right,# aspect=10,
        limits=(-1.2 * xmax, 1.2 * xmax + 2 * res, nothing, nothing), ylabel="Utilization [ratio]",
    )

    axs_moment_ratio = Axis(g[3, 1], yaxisposition=:right, #aspect=10,
        limits=(-1.2 * xmax, 1.2 * xmax + 2 * res, nothing, nothing), ylabel="Utilization [ratio]",
    )
    axs_shear_ratio = Axis(g[4, 1], yaxisposition=:right, #aspect=10,
        limits=(-1.2 * xmax, 1.2 * xmax + 2 * res, nothing, nothing), ylabel="Utilization [ratio]",
    )

    hidespines!(axs_axial_ratio)
    hidespines!(axs_moment_ratio)
    hidespines!(axs_shear_ratio)
    hideydecorations!(axs_axial_ratio, ticks=false, label=false, ticklabels=false)
    hideydecorations!(axs_moment_ratio, ticks=false, label=false, ticklabels=false)
    hideydecorations!(axs_shear_ratio, ticks=false, label=false, ticklabels=false)
    ypos = 0
    points = []
    if type[1] == 3
        for i in 1:n
            points = [x_range[i] - 250, -L, 500, L * 1.5]
            poly!(axs_design, Rect(points...), color=set_fc′[i],
                colorrange=(0, 80), colormap=cgrad(:Greys_3),
                strokecolor=:black, strokewidth=2)
            text!(axs_design, points[1], 100, text=string(set_fc′[i]))
            text!(axs_design, points[1], 100 + 75, text=string(set_dosage[i]))

        end
        section_pts = make_Y_layup_section(L, t, Lc)
        for s in section_pts.solids
            poly!(axs_design, s.points .+ [x_n + res + 500; 0], color=:grey)
        end
        ypos = 150
    elseif type[1] == 2
        for i in 1:n
            points = [x_range[i] - 250, -0.5 * L, 500, L]
            poly!(axs_design, Rect(points...), color=set_fc′[i],
                colorrange=(0, 80), colormap=cgrad(:Greys_3),
                strokecolor=:black, strokewidth=2)
            text!(axs_design, points[1], 100, text=string(set_fc′[i]))
            text!(axs_design, points[1], 100 + 75, text=string(set_dosage[i]))

        end
        #plot cross section on the side of the beam plot
        section_pts = make_X2_layup_section(L, t, Lc)
        for s in section_pts.solids
            poly!(axs_design, s.points .+ [x_n + res + 500; 0], color=:grey)

        end
        ypos = -250
    elseif type[1] == 4
        for i in 1:n
            points = [x_range[i] - 250, -L, 500, 2 * L]
            poly!(axs_design, Rect(points...), color=set_fc′[i],
                colorrange=(0, 80), colormap=cgrad(:Greys_3),
                strokecolor=:black, strokewidth=2)
            text!(axs_design, points[1], 100, text=string(set_fc′[i]))
            text!(axs_design, points[1], 100 + 75, text=string(set_dosage[i]))

        end
        #plot cross section on the side of the beam plot
        section_pts = make_X4_layup_section(L, t, Lc)
        for s in section_pts.solids
            poly!(axs_design, s.points .+ [x_n + res + 500; 0], color=:grey)
        end

    else
        println("Invalid Type")
    end

    Colorbar(f1[1, 2], ticks=[40, 60, 80], colorrange=(40, 100), colormap=cgrad(:Greys_3))


    tendon = poly!(axs_design, tendon_points, color=:skyblue, alpha=0.5, transparent=true)
    tendon_pts = scatter!(axs_design, tendon_points)
    tendon = lines!(axs_design, x_range, -tendon_profile, linewidth=5)
    tendon_constraints = lines!(axs_design, x_range,  -set_max_dps, color = :red, linewidth = 2)

    @show tendon_points
    for i in 1:n
        text!(axs_design, tendon_points[1, i+n], tendon_points[2, i+n], text=string(round(set_as[i], digits=2)))
        # string(round(set_axial_force[1], digits = 3))
    end


    #add text description
    txt_force = round(maximum(set_fps .* set_as) / 1000, digits=1) #kN.
    if type[1] == 2.0
        ypos1 = -200 
        ypos2 = -300
        ypos3 = -300
        xpos1 = -1000
        xpos2 = -1000
        xpos3 = 0
    else
        ypos = -400
        xpos = -xmax - 500

        ypos1= ypos; ypos2 = ypos; ypos3 = ypos; 
        xpos1 = xpos; xpos2 = xpos+1500; xpos3 = 0
    end

    text!(axs_design, xpos1, ypos1, text="Tendon load: " * string(txt_force) * " kN")
    # text!(axs_design,-xmax+700, ypos, text = "Apply test load "*string(set_P[1])* " kN")
    text!(axs_design, xpos2, ypos2, text="Axial post tensioned force " * string(round(set_axial_force[1], digits=1)) * " kN")
    text!(axs_design, xpos3, ypos3, text = "Carbon:"*string(round(sum(set_embodied_carbon), digits = 2)))


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


    axial_utilization = [axial_demand[i] / axial_capacity[i] for i in eachindex(axial_demand)]
    moment_utilization = [moment_demand[i] / moment_capacity[i] for i in eachindex(moment_demand)]
    shear_utilization = [shear_demand[i] / shear_capacity[i] for i in eachindex(shear_demand)]

    # moment_utilization =
    # shear_utilizatishear
    scatter!(axs_axial_ratio, x_range, axial_utilization, color=:red, marker=:utriangle)
    scatter!(axs_moment_ratio, x_range, moment_utilization, color=:blue, marker=:utriangle)
    scatter!(axs_shear_ratio, x_range, shear_utilization, color=:green, marker=:utriangle)



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