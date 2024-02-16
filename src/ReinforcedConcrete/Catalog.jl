# module Catalog
#design catalog
using DataFrames
using Makie, GLMakie
include("Definitions.jl")
include("Rebars/Rebars.jl")
include("RcCapacities.jl")

"""
β1 from ACI 318-19
"""
function get_β1(fc′::Real)
    return clamp(0.85 - 0.05 * (fc′ - 28) / 7, 0.65, 0.85)
end


"""
find the boundary (ρ_min, ρ_max) of the reinforcement ratio
ρmin is based on ACI318M-19 9.6.1.2
Input\\
fc′: concrete compressive strength [MPa]\\
Input [Optional]\\
fy = 420 steel yield strength [MPa]\

Outputs\\
ρ_min: minimum reinforcement ratio\\
ρ_max: maximum reinforcement ratio\\
"""
function get_ρ_bounds(fc′::Real; 
    fy::Float64 = 420.0)
    

    ρ_min1 =  0.25 * sqrt(fc′) / fy
    ρ_min2 =  1.4/fy
    ρ_min = maximum([ρ_min1, ρ_min2])

    β1 = get_β1(fc′)
    ρ_max = 3/5*0.85*β1*fc′/fy

    return ρ_min, ρ_max
end

"""
Generate the design space
note: Output as a dictionary might be better in term of performance.
"""
function get_design_space()
    set_fc′ = 25.0:0.5:55.0
    set_d   = 100.0:25:500.0
    set_bd_ratio = 0.5:0.05:1.0
    set_mu = 1e6(10:10:500) #[Nmm]
    # widths = 200.:100.:500.
    # heights = 200.:100.:500.
    # set_bd_ratio = 0.5:0.05:1 #b/d ratio b = bd_ratio*d
    # rebars = bar_combinations  #get from Hazel's work
    # rebars = bar_combinations  #get from azel's work
    

    # design_space = Dict()
    # for i1 in set_fc′
    #     for i2 in set_d
    #         for i3 in set_bd_ratio
    #             push!(design_space, (i1,i2,i3))
    #         end
    println("Design Space of ", 100*prod(length.([set_fc′, set_d, set_bd_ratio])) , " points")
    return set_fc′, set_d, set_bd_ratio, set_mu
end

function get_catalog()
    #some constants
    fy = 420.0
    covering = 40. #ACI318M-19 Table 20.5.1.3.1, Not exposed to weather or in contact with ground
     
    set_fc′, set_d, set_bd_ratio, set_mu = get_design_space()
    #results placeholder
    Cs = Vector{ConcreteSection}()
    Ps = Vector{Float64}() #Compressive Capacity
    Ms = Vector{Float64}() #Moment Capacity
    fc′s = Vector{Float64}()
    Areas = Vector{Float64}()
    ρs = Vector{Float64}()
    GWPs = Vector{Float64}()
    Section_IDs = Vector{Int64}()
    depths = Vector{Float64}()
    bd_ratios = Vector{Float64}() #b/d ratio

    #for different combination of sections and IDs.
    section_ID = 0
    count = 0 
    
    for bd_ratio in set_bd_ratio
        for d in set_d
            
            b = bd_ratio*d

            #reate an AsapSection section
            #rectangular section
            # p1.....p2
            # .      .
            # .      .
            # .      .
            # p4.....p3
            
            #for now w = b , d = h, normally w = b + 2*cover
            w = b
            h = d
            p1 = [0. , 0.]
            p2 = [w , 0.]
            p3 = [w, -h]
            p4 = [0., -h]
            pts = [p1,p2,p3,p4]
            section = SolidSection(pts)
            section_ID = section_ID + 1
            for fc′ in set_fc′
            # for mu in set_mu
                #get the bounds of ρ
                ρ_min , ρ_max = get_ρ_bounds(fc′)
                n_ρ = 100
                ρ_step = (ρ_max - ρ_min)/(n_ρ-1)
            
                #now, we will sample a 100 samples over the range of ρ
                for ρ in ρ_min:ρ_step:ρ_max
                    #calculate the As
                    as = ρ*b*d
                    #rectangular rebar (dummy)
                    rebars = RebarSection([as], [fy], [b/2], [-h + as/b], [0.])

                     
                    #Create a concrete section here.
                    # @show fc′ = find_fc′(mu, section, rebars)

                    c = ConcreteSection(fc′, section, rebars)
                    #Calculate the capacity.
                    P = find_Pu(c)
                    M = find_Mu(c)
                    if M < 0
                        continue
                    end
                    count = count +1
                    # V = find_Vn(c)
                    push!(Cs,c)
                    push!(Ps,P)
                    push!(Ms,M)
                    push!(GWPs, c.gwp)
                    push!(Section_IDs, section_ID)

                    #will get rid of this later
                    push!(fc′s, fc′)
                    push!(Areas, c.geometry.area)
                    push!(ρs, ρ)
                    push!(depths , d)
                    push!(bd_ratios, bd_ratio)
                end
            end
        end
    end

    # after got all of the catalog, compared them 
    catalog = DataFrame(id = 1:count, Section = Cs,fc′ = fc′s, Area = Areas, Depth = depths,bd_ratio = bd_ratios, ρ = ρs, Pu=Ps, Mu= Ms, Gwp = GWPs, Section_ID = Section_IDs)
    println("Done catalog with ", count, " points")
    return catalog
end

# end
catalog = get_catalog();
gwp_concrete = getfield.(catalog[!,:Section], :gwp_concrete)
gwp_rebars = getfield.(catalog[!,:Section], :gwp_rebars)
ratio_c = gwp_concrete ./ (gwp_rebars .+ gwp_concrete)
ratio_r = gwp_rebars ./ (gwp_rebars .+ gwp_concrete)

catalog[!, :ratio_c] = ratio_c
catalog[!, :ratio_r] = ratio_r
sorted_catalog = sort(catalog, [:Mu,:Gwp],rev = true)
pareto = sorted_catalog[1:1,:]

foreach(row -> row.Gwp < pareto.Gwp[end] && push!(pareto, row), eachrow(sorted_catalog));






##############
figure1 = Figure(resolution = (2000,1000)) 
ax1 = Axis(figure1[1,1], xlabel = "Moment Capacity [kN.m]", ylabel = "GWP [kgCO2e]")
s1 = scatter!(ax1, catalog.Mu/1e6, catalog.Gwp, markersize = 2,color = catalog.fc′,  label = "Catalog", strokewidth = 0)
Colorbar(figure1[2,1], s1, label = "fc′ [MPa]", labelrotation =0,vertical = false)

ax2 = Axis(figure1[1,2], xlabel = "Moment Capacity [kN.m]", ylabel = "GWP [kgCO2e]")
s1_2 = scatter!(ax2, catalog.Mu/1e6, catalog.Gwp, markersize = 2,color = :grey,  label = "Catalog", strokewidth = 0)
s2 = scatter!(ax2, pareto.Mu/1e6, pareto.Gwp, markersize = 10,color = pareto.fc′,  label = "Pareto", strokewidth = 0)
vlines!(ax2, 500, color = :red, label = "Pareto")
Colorbar(figure1[2,2], s2, label = "fc′ [MPa]", labelrotation =0,vertical = false,)





figure2 = Figure(resolution = (3000,1000))
# pareto_less_than_500 = pareto[pareto.Mu .< 500e6,:]
ax3 = Axis(figure2[1,1], xlabel = "Moment Capacity [kN.m]", ylabel = "GWP [kgCO2e]")
s2_3 = scatter!(ax3, pareto.Mu/1e6, pareto.Gwp, markersize = 10,color = pareto.fc′,  label = "Pareto", strokewidth = 0)
Colorbar(figure2[2,1], s2_3, label = "fc′ [MPa]", labelrotation =0,vertical = false,)

ax4 = Axis(figure2[1,2], xlabel = "Moment Capacity [kN.m]", ylabel = "GWP [kgCO2e]")
# s1_2 = scatter!(ax2, catalog.Mu/1e6, catalog.Gwp, markersize = 2,color = :grey,  label = "Catalog", strokewidth = 0)
s4 = scatter!(ax4, pareto.Mu/1e6, pareto.Gwp, markersize = 10,color = pareto.Depth,  label = "Pareto", strokewidth = 0)
# vlines!(ax, 500, color = :red, label = "Pareto")
Colorbar(figure2[2,2], s4, label = "Depth [mm]", labelrotation =0,vertical = false,)

#same thing but bd_ratio
ax5 = Axis(figure2[1,3], xlabel = "Moment Capacity [kN.m]", ylabel = "GWP [kgCO2e]")
# s1_2 = scatter!(ax2, catalog.Mu/1e6, catalog.Gwp, markersize = 2,color = :grey,  label = "Catalog", strokewidth = 0)
s5 = scatter!(ax5, pareto.Mu/1e6, pareto.Gwp, markersize = 10,color = pareto.bd_ratio,  label = "Pareto", strokewidth = 0)
# vlines!(ax, 500, color = :red, label = "Pareto")
Colorbar(figure2[2,3], s5, label = "bd_ratio", labelrotation =0,vertical = false,)

# save("limit at 500mm.png", figure2)
# save("limit at 1000mm.png",figure2)

# figure1

#Ratio plot 

figure3 = Figure(resolution = (3000,1000))
ax2 = Axis(figure3[1,1], xlabel = "Moment Capacity [kN.m]", ylabel = "c_ratio", limits = ( nothing, nothing,0,1))
ax1 = Axis(figure3[1,2], xlabel = "Moment Capacity [kN.m]", ylabel = "r_ratio",limits = ( nothing, nothing,0,1))
s_ratio = scatter!(ax1, pareto.Mu/1e6, pareto.ratio_c, color= :red)
s_r_ratio = scatter!(ax2, pareto.Mu/1e6, pareto.ratio_r, color = :blue)

save("Ratio.png", figure3)




"""
Visualize the design space
"""
function VizCatalog_section(catalog)
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

    #get for each Mu, get minimum gwp designs.

   
    figure1 = Figure(resolution = (600,600), backgroundcolor = :grey)
    ax1 = Axis(figure1[1,1], title = "Section Plot"
    ,xlabel = "Mu [kNm]", ylabel = "gwp [kgCO2e/kg.m]",
    limits = (0,Mu_max+10,0,1.1*gwp_max))
    slider1 = Slider(figure1[1,2], range = 1:maximum(catalog[!,:Section_ID]), startvalue =1,horizontal = false)
    
    x1 = Observable(catalog[!, :Mu]/1e6)
    y1 = Observable(catalog[!, :Gwp])
    c1 = Observable(catalog[!, :fc′])
    lift(slider1.value) do n
        new_cat = filter(:Section_ID => x-> x == n, catalog)
        # x[] = Point2f.(vcat.(new_cat[!,:area], new_cat[!,:gwp]))
        x1.val = new_cat[!,:Mu]/1e6
        y1[] = new_cat[!,:Gwp]
        # c1[] = new_cat[!,:fc′]
        c1[] = new_cat[!,:ρ]

        # @show x1.val[end]

        # z[] = new_cat[!,:fc′]
        # title_name[] = string(n)
        println(size(new_cat))
        return new_cat[!, :Mu], new_cat[!, :Gwp]
    end
    s0 = scatter!(ax1, catalog[!,:Mu]/1e6,catalog[!,:Gwp], color = :red, strokewidth=0, alpha = 0.1)
    s1 = scatter!(ax1, x1,y1, color = c1, strokewidth=0,alpha = 0.5)# colorrange = (fc′_min,fc′_max))
    Colorbar(figure1[2,1], s1, label = "fc′ [MPa]", labelrotation =0, vertical = false)

    # pairplot(catalog[!, [:Gwp]], catalog[!, [:fc′, :Area,:Mu, :Pu, :ρ]])
    return figure1
end
f_this = VizCatalog_section(catalog)