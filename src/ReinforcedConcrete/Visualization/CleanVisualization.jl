using Makie, GLMakie
using kjlMakie
set_theme!(kjl_light)
# using PairPlots
# using PlotlyJS

# figure1 = Figure(resolution = (1920, 1000))
# pairplot(catalog[!, [:Gwp, :Mu]], catalog[!, [:fc′, :Area,:Mu, :Pu, :ρ]] )
"""
Visualize the design space
"""
function VisCatalog(catalog)
    #set boundaries for plot
    Mu_min = minimum(catalog[!,:Mu])/1e6
    Mu_max = maximum(catalog[!,:Mu])/1e6
    fc′_min = minimum(catalog[!, :fc′])
    fc′_max = maximum(catalog[!, :fc′])
    area_min = minimum(catalog[!, :Area])
    area_max = maximum(catalog[!, :Area])
    ρ_min = minimum(catalog[!, :ρ])
    ρ_max = maximum(catalog[!, :ρ])
    gwp_min = minimum(catalog[!,:Gwp])
    gwp_max = maximum(catalog[!,:Gwp])


    # pairplot(catalog[!, [:Gwp]], catalog[!, [:fc′, :Area,:Mu, :Pu, :ρ]])
    figure1 = Figure(resolution = (1920, 1500),backgroundcolor = :grey)

    # gwp vs fc' (best)
    ax1 = Axis(figure1[1,1],
        xlabel = "fc′ [MPa]", ylabel = "GWP kgCO2e/kg", title = "Overall Plot",
        limits = (fc′_min-5,fc′_max+5,0,1.1*gwp_max))
    # fc′s = getfield.(catalog[!,:Section], :fc′)
    s1 = scatter!(ax1, catalog[!, :fc′], markersize = 5, 
    catalog[!, :Gwp],
    colormap =:amp, color = catalog[!, :Mu]/1e6,
    strokewidth = 0 )
    Colorbar(figure1[2,1], s1, label = "Mu [kNm]", labelrotation =0,vertical = false)

    ax2 = Axis(figure1[1,2],
    xlabel = "Mu [kNm]", ylabel = "GWP [kgCO2e/kg]",
    limits = (0,Mu_max+10,0,1.1*gwp_max))
    fc′s = getfield.(catalog[!,:Section], :fc′)
    s2 = scatter!(ax2, catalog[!, :Mu]/1e6, catalog[!, :Gwp], color = fc′s, markersize = 5 , strokewidth = 0)
    Colorbar(figure1[2,2], s2, label = "fc′ [MPa]", labelrotation =0, vertical = false)

    ax3 = Axis(figure1[1,3],
    xlabel = "SteelArea [mm2]", ylabel = "GWP [kgCO2e/kg]",
    limits = (0,area_max+10,0,1.1*gwp_max))
    areas = getfield.(getfield.(catalog[!,:Section], :geometry), :area)
    s3 = scatter!(ax3, areas, catalog[!, :Gwp], color = fc′s, markersize = 5 , strokewidth = 0)
    Colorbar(figure1[2,3], s3, label = "fc′ [MPa]", labelrotation =0, vertical = false)
    return figure1
end
f1 = VizCatalog(catalog)
f11 = VizCatalog(pareto)
# save("Monday_catalog1_1.png", f1)
#########################################
"""
Visualize the design space
"""
function VizCatalog_ratio(catalog)
    #set boundaries for plot
    Mu_min = minimum(catalog[!,:Mu])/1e6
    Mu_max = maximum(catalog[!,:Mu])/1e6
    fc′_min = minimum(catalog[!, :fc′])
    fc′_max = maximum(catalog[!, :fc′])
    area_min = minimum(catalog[!, :Area])
    area_max = maximum(catalog[!, :Area])
    ρ_min = minimum(catalog[!, :ρ])
    ρ_max = maximum(catalog[!, :ρ])
    gwp_min = minimum(catalog[!,:Gwp])
    gwp_max = maximum(catalog[!,:Gwp])

    gwp_concrete = getfield.(catalog[!,:Section], :gwp_concrete)
    gwp_rebars = getfield.(catalog[!,:Section], :gwp_rebars)
    ratio = gwp_concrete ./ (gwp_rebars .+ gwp_concrete)
    
    figure1 = Figure(resolution = (1920, 1000),backgroundcolor = :grey)

    ax1 = Axis(figure1[1,1],title = "Plot Cocncrete Ratio",
    xlabel = "Mu [kNm]", ylabel = "GWP [kgCO2e/kg]",
    limits = (0,Mu_max+10,0,1.1*gwp_max))
    s1 = scatter!(ax1, catalog[!, :Mu]/1e6, catalog[!, :Gwp], color = ratio, markersize = 10 , strokewidth = 0)
    Colorbar(figure1[2,1], s1, label = "Concrete contribution", labelrotation =0, vertical = false)

    ax2 = Axis(figure1[1,2],
    xlabel = "Mu [kNm]", ylabel = "Concrete contribution", #zlabel = "Concrete Contribution [ratio]",
    limits = (0,Mu_max+10, 0,1.1),)
    s2 = scatter!(ax2, catalog[!, :Mu]/1e6,ratio, color = ratio, markersize = 1 , strokewidth = 0)

    ax3 = Axis(figure1[1,3],
    xlabel = "Mu [kNm]", ylabel = "Concrete contribution", #zlabel = "Concrete Contribution [ratio]",
    limits = (0,Mu_max+10, 0,1.1),)
    s2 = scatter!(ax2, catalog[!, :Mu]/1e6,ratio, color = ratio, markersize = 1 , strokewidth = 0)



    return figure1
end

ratio_plot = VizCatalog_ratio(catalog)
ratio_plot2 = VizCatalog_ratio(pareto)


function VizCatalog_min(catalog)
    Mu_max = maximum(catalog[!,:Mu])/1e6
    gwp_max = maximum(catalog[!,:Gwp])

    sorted_catalog = sort(catalog, [:Mu,:Gwp],rev = true)
    paretoo = sorted_catalog[1:1,:]

    foreach(row -> row.Gwp < paretoo.Gwp[end] && push!(paretoo, row), eachrow(sorted_catalog));
    @show size(paretoo)
    @show paretoo[1,:]

    figure1 = Figure(resolution = (1920, 1000),backgroundcolor = :grey)
    ax1 = Axis(figure1[1,1],title = "Minimum plot",
    xlabel = "Mu [kNm]", ylabel = "GWP [kgCO2e/kg]",
    limits = (0,Mu_max+10,0,1.1*gwp_max))
    s1 = scatter!(ax1,
    paretoo[!, :Mu]/1e6 , 
    paretoo[!, :Gwp],
    color = paretoo[!, :Section_ID],
    strokewidth = 0 )
    Colorbar(figure1[2,1], s1, label = "fc′ [MPa]", labelrotation =0,vertical = false)

    return figure1, paretoo
end
f13d, pareto = VizCatalog_min(catalog)
# save("catalog13D.png", f13d)


save

f3 = VizCatalog_with_lines(catalog)





function VizCatalog_by_section(catalog)
    figure1 = Figure(resolution = (1920,1000))
    Section_IDs = unique(catalog[!,:Section_ID])
    ns = length(Section_IDs)
    ax1 = Axis(figure1[1,1],
        xlabel = "Mu [kNm]", ylabel = "GWP kgCO2e/kg")

    ax2 = Axis(figure1[1,2],
    xlabel = "Mu [kNm]", ylabel = "GWP kgCO2e/kg")
    fc′s = getfield.(catalog[!,:Section], :fc′)
    s1 = scatter!(ax1, catalog[!, :Mu]/1e6, catalog[!, :Gwp], color = fc′s, strokewidth = 0 )
    s2 = scatter!(ax2, catalog[!, :Mu]/1e6, catalog[!, :Gwp], color = catalog[!, :Area], #.+ fc′s./maximum(fc′s)/2,
    colormap = cgrad(:Spectral,ns, categorical = true), strokewidth = 0
    # opacity =  fc′s./maximum(fc′s)
    )
    Colorbar(figure1[2,1], s1, label = "fc′", vertical = false)
    Colorbar(figure1[2,2], s2, label = "Section", vertical = false)

    Label(figure1[1, 1:2, Top()], "Rc Beam Catalog", valign = :bottom,
    font = :bold, fontsize = 30,
    padding = (0, 0, 5, 0))
    # # #loop each section ID and plot it's lines
    # Section_IDs = unique(catalog[!,:Section_ID])
    # for i in eachindex(Section_IDs)
    #     section_ID = Section_IDs[i]
    #     selected = filter(:Section_ID => ==(section_ID), sort(catalog, [:Gwp, :Mu]))

    #     lines!(ax2, selected[:,:Mu]/1e6, selected[:, :Gwp])
    # # Colorbar(figure1[1,3], s2, label = "fc′", labelrotation =0)
    # end
    return figure1
end


f2 = VizCatalog_by_section(catalog)
# save("bySection.png", f2)



function VizCatalog_by_mu(catalog)
    figure1 = Figure(resolution = (1920,1000))

    #use menu to select mu range. 
    #a slider bar to select mu range, and show that on the plot
    Section_IDs = unique(catalog[!,:Section_ID])
    ns = length(Section_IDs)
    ax1 = Axis(figure1[1,1],
        xlabel = "fc′ [MPa]", ylabel = "GWP [kgCO2e/kg/m]")

    mumin = minimum(catalog[!,:Mu])/1e6
    mumax = maximum(catalog[!,:Mu])/1e6
    slider1 = Slider(figure1[1,2], range = mumin:mumax, startvalue = mumin,horizontal = false)
    slider2 = Slider(figure1[1,3], range = mumin:mumax, startvalue = mumax,horizontal = false)

    point = lift(slider1.value) do n
        Point2f(n,0)
    end

    label = lift(slider1.value, slider2.value) do s1, s2
        s1 = round(s1,digits =2)
        s2 = round(s2,digits =2)
        return "Range = $s1 to $s2"
    end
    Label(figure1[1,4], label, tellheight = false, fontsize = 30)

    minval = lift(slider1.value) do n
        return n*1e6
    end
    maxval = lift(slider2.value) do n
        return n*1e6
    end
    # minval = Observable{Any}(0.0)
    # maxval = Observable{Any}(0.0)
    # minval[] = slider1.value
    # maxval[] = slider2.value

    #new dataframe 
    filtered_df = filter(:Mu => x-> x > to_value(minval) && x < to_value(maxval), catalog)
    df_obs = Observable(catalog)
    data = @lift(Point2f.($df_obs[:, $col_1], $df_obs[:, $col_2]))

    # text!(0.5, 0.5, text = String(MuVal))
    fc′s = Observable(getfield.(filtered_df[!,:Section], :fc′))
    s1 = scatter!(ax1, fc′s, Observable(filtered_df[!, :Gwp]), color = :blue)
    s2 = scatter!(ax1, point, color = :red)
    # ax2 = Axis(figure1[1,2],
    # xlabel = "Mu [kNm]", ylabel = "GWP kgCO2e/kg")

    # s1 = scatter!(ax1, filtered_df[!, :Mu]/1e6, filtered_df[!, :Gwp], color = :grey)
    # s2 = scatter!(ax2, catalog[!, :Mu]/1e6, catalog[!, :Gwp], color = catalog[!, :Section_ID],
    # colormap = cgrad(:Spectral,ns, categorical = true))
    # Colorbar(figure1[2,1], s1, label = "fc′", vertical = false)
    # Colorbar(figure1[2,2], s2, label = "Section", vertical = false)

    # Label(figure1[1, 1:2, Top()], "Rc Beam Catalog", valign = :bottom,
    # font = :bold, fontsize = 30,
    # padding = (0, 0, 5, 0))
    # # #loop each section ID and plot it's lines
    # Section_IDs = unique(catalog[!,:Section_ID])
    # for i in eachindex(Section_IDs)
    #     section_ID = Section_IDs[i]
    #     selected = filter(:Section_ID => ==(section_ID), sort(catalog, [:Gwp, :Mu]))

    #     lines!(ax2, selected[:,:Mu]/1e6, selected[:, :Gwp])
    # # Colorbar(figure1[1,3], s2, label = "fc′", labelrotation =0)
    # end
    return figure1
end

function VizCatalog_by_mu(catalog, minval, maxval)
    figure1 = Figure(resolution = (1920,1000))

    #use menu to select mu range. 
    #a slider bar to select mu range, and show that on the plot
    Section_IDs = unique(catalog[!,:Section_ID])
    ns = length(Section_IDs)
    ax1 = Axis(figure1[1,1],
        xlabel = "fc′ [MPa]", ylabel = "GWP [kgCO2e/kg/m]",
        limits = (20,60,0,120),
        xlabelsize = 30, ylabelsize = 30,
        xminorticksize = 30, title ="Range = $minval to $maxval kNm", titlesize = 40
        )

    # mumin = minimum(catalog[!,:Mu])/1e6
    # mumax = maximum(catalog[!,:Mu])/1e6
    # slider1 = Slider(figure1[1,2], range = mumin:mumax, startvalue = mumin,horizontal = false)
    # slider2 = Slider(figure1[1,3], range = mumin:mumax, startvalue = mumax,horizontal = false)

    # point = lift(slider1.value) do n
    #     Point2f(n,0)
    # end

    # label = lift(slider1.value, slider2.value) do s1, s2
    #     s1 = round(s1,digits =2)
    #     s2 = round(s2,digits =2)
    #     return "Range = $s1 to $s2"
    # end
    # Label(figure1[1,2], label, tellheight = false, fontsize = 30)
    # minval = lift(slider1.value) do n
    #     return n*1e6
    # end
    # maxval = lift(slider2.value) do n
    #     return n*1e6
    # end
    # minval = Observable{Any}(0.0)
    # maxval = Observable{Any}(0.0)
    # minval[] = slider1.value
    # maxval[] = slider2.value

    #new dataframe 
    filtered_df = filter(:Mu => x-> x > minval*1e6 && x < maxval*1e6, catalog)

    # text!(0.5, 0.5, text = String(MuVal))
    fc′s = getfield.(filtered_df[!,:Section], :fc′)
    s1 = scatter!(ax1, fc′s, filtered_df[!, :Gwp], color = filtered_df[!, :Mu],
    markersize = 20,
    inspector_label = (self, i, p) -> string(filtered_df[i, :Section_ID])
    )
    Colorbar(figure1[2,1], s1, label = "Mu", labelrotation =0, vertical = false)

    # s2 = scatter!(ax1, point, color = :red)
    # ax2 = Axis(figure1[1,2],
    # xlabel = "Mu [kNm]", ylabel = "GWP kgCO2e/kg")

    # s1 = scatter!(ax1, filtered_df[!, :Mu]/1e6, filtered_df[!, :Gwp], color = :grey)
    # s2 = scatter!(ax2, catalog[!, :Mu]/1e6, catalog[!, :Gwp], color = catalog[!, :Section_ID],
    # colormap = cgrad(:Spectral,ns, categorical = true))
    # Colorbar(figure1[2,1], s1, label = "fc′", vertical = false)
    # Colorbar(figure1[2,2], s2, label = "Section", vertical = false)

    # Label(figure1[1, 1:2, Top()], "Rc Beam Catalog", valign = :bottom,
    # font = :bold, fontsize = 30,
    # padding = (0, 0, 5, 0))
    # # #loop each section ID and plot it's lines
    # Section_IDs = unique(catalog[!,:Section_ID])
    # for i in eachindex(Section_IDs)
    #     section_ID = Section_IDs[i]
    #     selected = filter(:Section_ID => ==(section_ID), sort(catalog, [:Gwp, :Mu]))

    #     lines!(ax2, selected[:,:Mu]/1e6, selected[:, :Gwp])
    # # Colorbar(figure1[1,3], s2, label = "fc′", labelrotation =0)
    # end
    #hover over the point, and get ID
    return figure1
end

mu = 50
ϵ  = 0.25
f2 = VizCatalog_by_mu(catalog,mu, mu+ϵ)
DataInspector(f2, indicator_color = :red,enable_indicators = true)


function ParplotCatalog(df) 
    fc′s = getfield.(catalog[!,:Section], :fc′)
    mytrace = parcoords(;line = attr(color=df.Gwp)
                 ,dimensions = [ attr(
                    label = "fc′",
                    values = fc′s,
                 ),
                            attr(
                                label = "Concrete gross area",
                                values = getfield.(getfield.(catalog[!,:Section], :geometry),:area)
                                ),
                            # attr(label = "Section Depth",
                            #     values = 
                            #     ),
                            attr(label = "GWP",
                                 values = df.Gwp
                                ),
                            attr(label = "Mu",
                            values = df.Mu
                            ),
                            attr(label = "Pu",
                            values = df.Pu),
                            # attr(label = "Steel Area",
                            # values = )
                   ]);

    layout = Layout(title_text="Parallel Coordinates Plot"
                  , title_x=0.5 
                  , title_y=0
                   )

    myplot = PlotlyJS.plot(mytrace,layout)

    return myplot
end


# open("./parallel.html", "w") do io
#     PlotlyBase.to_html(io, ParplotCatalog(catalog).plot)
# end