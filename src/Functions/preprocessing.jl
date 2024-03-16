function preprocessing(demand_path::String)::DataFrame

    println("Plot saved at " * imagesavepath *"*")
    open(demand_path, "r") do f
        global demands = DataFrame(JSON.parse(f, dicttype=Dict{String,Any}))
        ns = size(demands)[1]
        demands[!, :idx] = 1:ns #add indices into the entry.
        println("There are $ns points")
        println("Demands were loaded from:\n", demand_path)
    end
    println(demands[1:10, :])

    #make the correct types
    demands[!, :e_idx]= Int.(demands[!,:e_idx])
    demands[!, :s_idx]= Int.(demands[!,:s_idx])
    demands[!, :idx]= Int.(demands[!,:idx])
    demands[!, :ec_max]= Float64.(demands[!,:ec_max])
    demands[!, :pu]= Float64.(demands[!,:pu])
    demands[!, :mu]= Float64.(demands[!,:mu])
    demands[!, :vu]= Float64.(demands[!,:vu])



    #Modify the section and element indices (if needed)
    println("__________")

    println("Before Modifying the indices")
    @show old_min_e_idx = minimum(demands[!, "e_idx"])
    @show old_min_s_idx = minimum(demands[!, "s_idx"])
    # ==================================================
    @show old_max_e_idx = maximum(demands[!, "e_idx"])
    @show old_max_s_idx = maximum(demands[!, "s_idx"])

    if old_min_e_idx == 0
        demands[!, "e_idx"] .+= 1
    else
        println("***Element indices are not modified")
    end

    if old_min_s_idx == 0
        demands[!, "s_idx"] .+= 1
    else
        println("***Section indices are not modified")
    end
    println("__________")

    println("After Modifying the Indices")
    @show new_min_e_idx = minimum(demands[!, "e_idx"])
    @show new_min_s_idx = minimum(demands[!, "s_idx"])
    println("__________")

    @show new_max_e_idx = maximum(demands[!, "e_idx"])
    @show new_max_s_idx = maximum(demands[!, "s_idx"])
    println("__________")

    """
Make the element indices unique
"""
    global e_idx = 1
    for i in 1:size(demands)[1]
        if i != 1
            demands[i-1, :e_idx] = e_idx
            if demands[i, :s_idx] < demands[i-1, :s_idx]
                global e_idx += 1
            end
        end
        if i == size(demands)[1]
            demands[i, :e_idx] = e_idx
        end
    end
    #add colors to the types of sections ()
    types = unique(demands[!, :type])
    color_mapping = Dict(types .=> 1:1:length(types))

    f_indices_check = Figure(size=(800, 500))
    ax_indices_check = Axis(f_indices_check[1, 1], title="Indices Check", titlesize=40,
        xlabel="Global Index", ylabel="Element Index", xticks=0:10:size(demands)[1]*1.2, yticks=1:1:30)

    scatter!(ax_indices_check, Int.(demands[!, "idx"]), Int.(demands[!, "e_idx"]), color=[color_mapping[t] for t in demands[!, :type]])
    f_indices_check

    save(imagesavepath * "f_indices_check.png", f_indices_check)

    f_indices_check_mod = Figure(size=(1500, 1500), backgroundcolor=:white)
    ax_indices_check_mod = Axis(f_indices_check_mod[1, 1], xlabel="Global Index", ylabel="Element Index", xticks=0:20:250, yticks=1:1:30,
        title="Modified element indices")
    ax_section_indices_check_mod = Axis(f_indices_check_mod[1, 2], xlabel="Global Index", ylabel="Element Index", xticks=0:20:250, yticks=1:1:30,
        title="Section indices")
    ax_dps = Axis(f_indices_check_mod[2, :], xlabel="Global Index", ylabel="dps", xticks=0:20:250,
        title="Section indices")

    ax_m = Axis(f_indices_check_mod[3, :], xlabel="Global Index", ylabel="Moment", xticks=0:20:250,
        title="Section indices")
    ax_p = Axis(f_indices_check_mod[4, :], xlabel="Global Index", ylabel="Axial", xticks=0:20:250,
        title="Section indices")
    ax_v = Axis(f_indices_check_mod[5, :], xlabel="Global Index", ylabel="Shear", xticks=0:20:250,
        title="Section indices")

    scatter!(ax_indices_check_mod, demands[!, "idx"], demands[!, "e_idx"], color=[color_mapping[t] for t in demands[!, :type]])
    scatter!(ax_section_indices_check_mod, demands[!, "idx"], demands[!, "s_idx"], color=[color_mapping[t] for t in demands[!, :type]])
    scatter!(ax_dps, demands[!, "idx"], demands[!, "ec_max"] .* 1000, color=[color_mapping[t] for t in demands[!, :type]])
    scatter!(ax_m, demands[!, "idx"], demands[!, "mu"] .* 1000, color=[color_mapping[t] for t in demands[!, :type]])
    scatter!(ax_p, demands[!, "idx"], demands[!, "pu"] .* 1000, color=[color_mapping[t] for t in demands[!, :type]])
    scatter!(ax_v, demands[!, "idx"], demands[!, "vu"] .* 1000, color=[color_mapping[t] for t in demands[!, :type]])

    f_indices_check_mod
    save(imagesavepath * "f_indices_check_mod.png", f_indices_check_mod)

    return demands
end


function plot_distribution(demands, catalog)
    """
Visualize the demand points
"""
f1 = Figure(size = (2000,1500), backgroundcolor = :white)
mu_max = 1.2*maximum(demands[!, :mu])
vu_max = 1.2*maximum(demands[!, :vu])

ax1 = Axis(f1[1,1:2], xlabel = "Moment [kNm]", ylabel = "Shear [kN]", limits = (-10,nothing,-10,nothing), title = "Demand vs Catalog Space") 
# ax2 = Axis(f1[1,2], xlabel = "Moment [kNm]", ylabel = "Shear [kN]", title= "Catalog ONLY space")
	
ax3 = Axis(f1[2,1], xlabel = "Shear [kN]" , 
# limits = (-10, vu_max, nothing, nothing), 
title = "Shear (Demands)")
ax4 = Axis(f1[3,1], xlabel = "Shear [kN]" , 
# limits = (-10, vu_max, nothing, nothing), 
title = "Shear (Catalog)")

ax5 = Axis(f1[2,2], xlabel = "Moment [kNm]", title ="Moment (Demands)")
ax6 = Axis(f1[3,2], xlabel = "Moment [kNm]", title ="Moment (Catalog)")


scatter!(ax1, catalog[!, :Mu], catalog[!,:Vu], color= "#7acfff", markersize = 5)
scatter!(ax1, demands[!, :mu], demands[!,:vu], color = "#911e98", label = demands[!, :type], marker = 'x')
# scatter!(ax2, catalog[!, :Mu], catalog[!, :Vu] , color = collect(catalog[!, :fc′]), markersize = 5)

hist!(ax3, demands[!, :vu], color = :red, bins = 20)
hist!(ax4, catalog[!, :Vu], color = :blue, bins = 20)
hist!(ax5, demands[!, :mu], color = :red, bins = 20)
hist!(ax6, catalog[!, :Mu], color = :blue, bins = 20)

elem_1 = [PolyElement(color = :red, strokecolor = :grey, strokewidth = 1)]
elem_2 = [PolyElement(color = :blue, strokecolor = :grey, strokewidth = 1)]
Legend(f1[2, 2],
    [elem_1, elem_2],
    ["Demand", "Catalog"],
    patchsize = (35, 35), rowgap = 10,
    halign = :right, valign = :top,
    tellheight = false,tellwidth = false)


f1
save(imagesavepath*"f_demand_catalog_distribution.png", f1)

end
