### A Pluto.jl notebook ###
# v0.19.38

using Markdown
using InteractiveUtils

# ╔═╡ 4f705574-ec7b-4243-ae03-1da701edac30
begin
    include("Functions/definition.jl");
    include("Functions/functions.jl");
    include("Functions/structuralelement.jl")
    include("Functions/generalfunctions.jl")
    
    include("Functions/get_Deflection.jl")
    include("Functions/interpolations.jl");
end

# ╔═╡ 8e24d686-a754-40e3-85fe-39acd8e94ae6
md"""
# PixelFrame engineering design demo"""

# ╔═╡ 0b2b7afa-f90b-4367-9798-1348128ddc5b
md"""
To do

actual length for absolute value of tendon profile
servicability from the model<Implemented but stuck in a loop"""

# ╔═╡ 0faeda67-c89e-49fe-bb4e-d5da4684b9b4
md"""
Import required packages"""

# ╔═╡ 8b6bfef7-cf3f-4ff0-8c55-f1ac29cf1525
# ╠═╡ skip_as_script = true
#=╠═╡
using Pkg
Pkg.activate("..");
Pkg.status()
Pkg.instantiate()

  ╠═╡ =#

# ╔═╡ d45ec2fe-84ea-4df6-a91a-d2e823249aad
md"""
Load functions associate with deflection calculation"""

# ╔═╡ 378a9fb6-d030-40c5-8e67-8a3e8b17697b
md"""
Check current directory"""

# ╔═╡ e2e244f6-ec5e-42e5-b3b6-e0fabad60630
pwd()

# ╔═╡ dc282072-4a34-485f-89f9-cc250ad93e6e
md"""
Load the precalculated catalog"""

# ╔═╡ fb79a319-57e8-4eec-b8c3-e392983fed7c
begin
    # catalog = CSV.read("Catalogs/test_catalog.csv", DataFrame); #load the pre-calc catalog
    catalog = CSV.read("Catalogs/FEB1_1_catalog_static.csv", DataFrame); #load the pre-calc catalog
    
    sort!(catalog, [:carbon, :fc′, :as, :dps])
    
    
    println("The catalog was sorted by ascending order from:\ncarbon -> fc′ -> as -> dps")
    println(catalog[1:20,:])
end

# ╔═╡ 4e048990-86f5-4bb1-86d3-7f04e5a81d8b
md"""
Load the demand file   """

# ╔═╡ 748c0454-6c73-4235-b957-92d55a29001e
md"""
The demand catalog was labeled with 0 index scheme, add 1 so they are compatible with Julia 1 index scheme."""

# ╔═╡ a15f7dd8-5711-46ad-a01a-8911cc17aa07
begin
    @show old_min_e_idx = minimum(demands[!, "e_idx"]);
    @show old_min_s_idx = minimum(demands[!, "s_idx"]);
    
    @show old_max_e_idx = maximum(demands[!, "e_idx"]);
    @show old_max_s_idx = maximum(demands[!, "s_idx"]);
    

end

# ╔═╡ b75a66bc-ae7e-446a-83e0-66cfd1d3fdd4
begin
    demands[!,"e_idx"] .+= 1 ;
    demands[!,"s_idx"] .+= 1 ;
end

# ╔═╡ 6572537d-7159-4bf2-9053-789327d7dad7
begin
    @show new_min_e_idx = minimum(demands[!, "e_idx"]);
    @show new_min_s_idx = minimum(demands[!, "s_idx"]);
    
    @show new_max_e_idx = maximum(demands[!, "e_idx"]);
    @show new_max_s_idx = maximum(demands[!, "s_idx"]);

end

# ╔═╡ 754f5262-6982-4779-991b-ed592ddc4358
begin
    f1 = Figure(size = (500,500))
    ax1 = Axis(f1[1,1], xlabel = "Moment [kNm]", ylabel = "Shear [kN]", limits = (0,3000,0,300)) 
    types = unique(demands[!, :type])
    @show mapping = Dict(types .=> 1:1:length(types))
    scatter!(ax1, demands[!, :mu], demands[!,:vu], color = [mapping[t] for t in demands[!,:type]], label = demands[!, :type])
    scatter!(ax1, catalog[!, :Mu], catalog[!, :Vu], marker = '.', alpha = 0.5, transparency = true, markersize = 5)
    f1
end

# ╔═╡ 3db397c5-dea4-4bf1-b431-e2af9f6ba53d
md"""
Fix redundant element indices due to Primary and Secondary labels"""

# ╔═╡ 17db0fd9-3260-46ed-a7b6-26e2a11225c4
begin
    global e_idx = 1 
    for i in 1:size(demands)[1]
        if i !=1
            demands[i-1, :e_idx] = e_idx
            if demands[i,:s_idx] < demands[i-1,:s_idx]
                global e_idx +=1 
            end
            
        end
        if i == size(demands)[1]
            demands[i, :e_idx] = e_idx
        end
    end
end

# ╔═╡ d5e7772a-c241-4e4a-8bc0-5d5ebefddb1f
md"""
Create a function to filter only feasible designs in a demand point"""

# ╔═╡ 1f7cfa53-4a8a-4d3c-932c-6c3d31a45d65
begin
    """
    filter feasible demands
    """
    function filter_demands(demands::DataFrame, catalog::DataFrame)::Dict{Int64, Vector{Int64}}
    
        ns = size(demands)[1]          #total number of sections
        ne = unique(demands[!, :e_idx]) #total number of elements
        nc = size(catalog, 1)           #total number of available choices
    
        all_feasible_sections = Dict{Int64,Vector{Int64}}() #map between demands and indices of the feasible section.
        demands[!, "total_results"] = zeros(Int64, size(demands)[1])
        #go through each section and filter the feasible designs from the catalog.
        for i = 1:ns
            en = demands[i, "e_idx"]
            sn = demands[i, "s_idx"]
            # push!(elements_to_sections[en], sn)
    
            pu = demands[i, "pu"]
            mu = demands[i, "mu"]
            vu = demands[i,"vu"]*0.9
            ec_max = demands[i, "ec_max"]
    
            global feasible_sections = filter([:Pu, :Mu, :Vu, :dps] => (x1, x2, x3, x4) ->
                    x1 > pu &&
                    x2 > mu &&
                    # x3 > vu && 
                    x4 <= ec_max * 1000,
                catalog
            )
    
            if size(feasible_sections)[1] == 0 #if the number of feasible results = 0
                println(feasible_sections[!, :ID])
                println("section $sn: element $en")
                all_feasible_sections[i] = [0]
                demands[i, "total_results"] = 0
                # println(outr)
            else
                all_feasible_sections[i] = feasible_sections[!, :ID]
                demands[i, "total_results"] = length(feasible_sections[!, :ID])
                # println(outr)
            end
        end
    
        return all_feasible_sections
    end

end

# ╔═╡ 70f74da6-d4b2-43aa-b3cf-be8aceff4220
begin
    """
    Find the optimum result for each element.
        
        !! not optimum yet.
    For the same element, will use the same fc′ steel size and post tensioning stress.
    
    """
    function find_optimum(all_feasible_sections::Dict{Int64, Vector{Int64}}, demands::DataFrame)
        total_ns = size(demands)[1] #get total number of section points.
        ne = unique(demands[!, :e_idx]) 
    
        #Map element indices to indices of sections indices on those elements.
        elements_to_sections = Dict(k => Int[] for k in unique(demands[!, :e_idx]))
        for i = 1:total_ns
            en = demands[i, :e_idx]
            # sn = demands[i, "s_idx"]
            push!(elements_to_sections[en], i)
        end
    
        #element index to list of designs
        elements_designs = Dict(k => [[]] for k in unique(demands[!, :e_idx])) 
        sections_to_designs = Dict(k => Vector{Float64}() for k in demands[!, :idx])
    
        #find fc′, as, and fpe that appear in all sections in an element.
        # for i in ProgressBar(ne) #loop each element
        for i in ne
            println("Element $i out of $(length(ne)) elements")
            sections = elements_to_sections[i] #sections associated with this element
    
            ns = length(sections)     # number of sections in this element 
    
            #start from the middle-ish section (n/2 or (n-1)/2)
            #note that section is in the form of 1,2,3,..., ns.
            mid = div(ns, 2)
            
            #get the feasible designs for the middle section
            feasible_idx = all_feasible_sections[sections[mid]]
    
            #catalog was already sorted, so I think we can leave this part, just filter, to save time.
            mid_catalog = sort(catalog[feasible_idx, :], [:carbon, :fc′, :dps])
            # mid_catalog = catalog[feasible_idx, :]
    
            #now, loop each design in the sub catalog, see if "as" and "fpe" are available in all sections.
            #if not, remove that design from the sub catalog.
            #if yes, keep it.
    
            #select each design, check if as and fpe exist for the the section
            global_d = 0 
            for d_idx in 1:size(catalog)[1]
                d = mid_catalog[d_idx, :]
                println(d)
    
                all_as  = true
                all_fpe = true
                serviceability_check = true
    
                for s in sections
                    #check if the design is available in that section.
                    #if not, remove it from the sub catalog.
                    #if yes, keep it.
    
                    #at this version, make sure that the fR1 and fR3 are also matched.
                    # if !(d[:fc′] ∈ catalog[all_feasible_sections[s], :fc′])
                    #     all_fc′ = false
                    #     break
                    # end
    
                    if !(d[:as] ∈ catalog[all_feasible_sections[s], :as])
                        all_as = false
                        break
                    end
                    if !(d[:fpe] ∈ catalog[all_feasible_sections[s], :fpe])
                        all_fpe = false
                        break
                    end
                end
            end
    
            if !all_fpe ||!all_as 
                println("Warning, can't find the solution for element $i")
            end
            
    
            #get the first one, they will appear in the entire thing anyway.
            this_fc′ = mid_catalog[global_d, :fc′]
    
            this_fpe = mid_catalog[global_d, :fpe]
            this_as = mid_catalog[global_d, :as]
    
            sections_designs = Vector{Vector}(undef, ns)
            for is in eachindex(elements_to_sections[i])
                #current section index
                s = elements_to_sections[i][is]
    
                feasible_idx = all_feasible_sections[s]
                
                # fc′_fpe_as(fc′::Float64, fpe::Float64, as::Float64) = fc′ == this_fc′ && fpe == this_fpe && as == this_as
                fpe_as(fpe::Float64, as::Float64) = fpe == this_fpe && as == this_as
    
                # this_catalog = filter([:fc′, :fpe, :as] => fc′_fpe_as, catalog[output_results[s], :])
                this_catalog = filter([:fpe, :as] => fpe_as, catalog[output_results[s], :])
    
                sort!(this_catalog, [:carbon, :dps])
    
                #get the first one, it's the lowest carbon
                select_ID = this_catalog[1, :ID]
                #find lowest e for this one.
                sections_designs[is] = collect(catalog[select_ID, :])
                # println(sections_designs[is])
            end
    
            #create a pixelframeelement and/or section here with the given parameters 
            L, t, Lc = [205.0 35.0 30.0]
            compoundsection =  make_Y_layup_section(L, t, Lc)
            pixelframeelement = PixelFrameElement() 
            Le = ns*500.0
    
            #could do a case where input only the variables -> the dependent variables come later.
            pixelframeelement.fc′ = this_fc′ # Concrete strength [MPa] ****Should update on the test day using cylinder test***
            # Ec = 4700.0*sqrt(fc′) # MPa  ACI fc-> Concrete modulus relationship [MPa]
            pixelframeelement.Ec = 58000.0 # MPa  from the cylinder test
            pixelframeelement.Eps = 70000.0 #Post tensioning steel modulus [MPa]
            pixelframeelement.fpy = 0.002 * pixelframeelement.Eps #MPa  
            #Safe load on the website https://www.engineeringtoolbox.com/wire-rope-strength-d_1518.html 
            # is ~ 150 MPa. Currently 140 MPa :)
    
            # PixelFrame section/element properties
            centroid_to_top = 100.0 #[mm] ~half of 205mm
            pixelframeelement.em = mid_catalog[global_d, :dps]# Eccentricity at the middle of the member [mm]
            pixelframeelement.es = 0.0 # Eccentricity at the support of the member   [mm]
            pixelframeelement.em0 = mid_catalog[global_d, :dps] # Initial eccentricity at the midspan        [mm]
    
            pixelframeelement.dps0 = centroid_to_top + pixelframeelement.em0 # Initial distance from the top to the point of application of the load [mm]
            pixelframeelement.Ls = Le/4 # Distance from support to the first load point [mm]
            pixelframeelement.Ld = Le/4 # Distance from support to the first deviator [mm]
            pixelframeelement.L = Le # Total length of the member [mm]
            # two 1/4" bars with 1200 lb capacity
            pixelframeelement.Aps = this_as # Total area of the post tensioned steel [mm2]
            #Pure concrete area = 18537.69 mm2
            #Transformed steel area = 347.96 mm2 
            pixelframeelement.Atr = compoundsection.area # Transformed area of the cross section [mm2] (= Concrete area if there is no embedded rebars)
            pixelframeelement.Itr = compoundsection.Ix #moment of inertia [mm4], no embedded steel, therefore, only from concrete.
            # pixelframeelement.Itr = 1.082e+8 #this number includes deviated steels.
    
            pixelframeelement.Zb = pixelframeelement.Itr/centroid_to_top # Elastic modulus of the concrete section from the centroid to extreme tension fiber [mm3]
            # If there are multiple materials, transformed section geometry is needed for Zb (and everything related to section area)
    
            #forces
            pixelframeelement.w = pixelframeelement.Atr / 10^9 * 2400.0 * 9.81 # Selfweight [N/mm]
            pixelframeelement.mg = pixelframeelement.w * pixelframeelement.L^2 / 8.0 # Moment due to selfweight [Nmm]
            pixelframeelement.fr = 0.7 * sqrt(pixelframeelement.fc′) # Concrete cracking strenght [MPa]
            pixelframeelement.r = sqrt(pixelframeelement.Itr / pixelframeelement.Atr) # Radius of gyration [mm]
            pixelframeelement.ps_force = pixelframeelment.Aps*this_fpe # Post tensioning force [N]
            pixelframeelement.Mdec = pixelframeelement.ps_force*pixelframeelement.em
            pixelframeelement.concrete_force = pixelframeelement.ps_force*cos(24.0*pi/180.0) # should use actual value 
            pixelframeelement.fpe = pixelframeelement.ps_force/pixelframeelement.Aps # Effective post tensioning stress [MPa] ***will input the one on the test day***
            pixelframeelement.ϵpe = pixelframeelement.fpe / pixelframeelement.Eps # Effective post tensioning strain [mm/mm]
            #find moment due to the applied force.
            pixelframeelement.ϵce = pixelframeelement.ps_force*pixelframeelement.em/pixelframeelement.Zb/pixelframeelement.Ec - pixelframeelement.concrete_force/pixelframeelement.Atr/pixelframeelement.Ec # effetive strain in the concrete [mm/mm]
            #for using test setup
            pixelframeelement.test = false
    
            #for the current stage of the model,
            #we will have 2 deviators, with interpolated 2/3 of the length of the element.
            #with this, <Need proved, but this should be more conservative???>. Less post tensioned than the actual beam.
            #Current PixelFrameSize
            
            load_m = 4*demands[mid, :mu]/Le*1000 #N
            load_v = demands[mid, :vu]*1000
    
            load = load_m > load_v ? load_m : load_v
            @show load 
            break
            dis_history, P =  get_Deflection(pixelframeelement, load, loadstep = 1000)
    
            δ = dis_history[end]
    
            elements_designs[i] = [sections_designs, δ, δ/(Le/240)]
    
        end
    
        for i in ne
            sections = elements_to_sections[i]
            for design in eachindex(elements_designs[i])
                sections_to_designs[sections[design]] = elements_designs[i][design]
            end
        end
    
        return elements_designs, elements_to_sections, sections_to_designs
    end
end

# ╔═╡ a1e47f88-ebd0-4773-88df-1ca6fc729e9c
md"""
Filter the catalog entries that meet the demands."""

# ╔═╡ 1c6d6600-88b6-4dfc-980e-deb6e18bb0c8
all_feasible_sections= filter_demands(demands,catalog)


# ╔═╡ f65446a9-ea2c-4443-9aff-354c61452bd5
md"""
Check number of available designs per section."""

# ╔═╡ ea1d4559-5f72-4c66-b2b3-62babc8cf708
# elements_designs, elements_to_sections  = find_optimum(output_results, demands)
elements_designs, elements_to_sections, sections_to_designs  = find_optimum(all_feasible_sections, demands)



# ╔═╡ af95e84f-c3c1-41c4-b439-e9a7ee6475d8
begin
    sections_to_designs
    open("Results/0102_section_to_designs.json","w") do f
        JSON.print(f, sections_to_designs)
    end
end

# ╔═╡ 54e76ba0-6e02-4f51-ba4a-a21838f6c495
@show sections_to_designs

# ╔═╡ 725672e0-7ae6-43da-90f8-98b3cc1911f2
md"""
Save the design result"""

# ╔═╡ b49ed097-8ef7-4187-bfd9-f7a4400c5e30
open("Results/designs_results_22_01.json","w") do f
    JSON.print(f, elements_designs)
end

# ╔═╡ 64a9cf05-2327-41c4-bc9f-19ff33e20049
md"""
## 3. Visualizing the results"""

# ╔═╡ 2ec7f3d1-9d2d-4ee9-9986-1218bbd22839
md"""
Load required packages"""

# ╔═╡ 8d6daa5d-ea86-45f1-bb0c-a12d13d34618
md"""
Load back the design results"""

# ╔═╡ d548ba9d-47ec-4f81-bb0e-94d9eb1843b9
begin
    designs = JSON.parsefile(joinpath(@__DIR__,"Results/designs_results_22_01.json"), dicttype = Dict{String,Vector{Vector{Float64}}});
    
    ne = length(designs)
    println("There are $ne elements.")

end

# ╔═╡ 8eed49de-fd5e-45ad-8bab-3a590bd1b7bd
designs["1"]

# ╔═╡ 5c26ac69-fd92-4b24-899e-890b860d9ade
designs["5"]

# ╔═╡ a7c95112-1348-46ef-b973-72ebcfb3276a
function plot_element(element_number::Int64, designs::Dict;L::Float64 = 250.0)


sections = elements_to_sections[element_number]
element_number = string(element_number)

L = 400
tendon_profile = L.*[i[3] for i in designs[element_number]]
axial_capacity  = [i[5] for i in designs[element_number]]
moment_capacity = [i[6] for i in designs[element_number]]
shear_capacity = [i[7] for i in designs[element_number]]


axial_demand  = [ demands[i,"pu"] for i in sections]
moment_demand = [ demands[i,"mu"] for i in sections]
shear_demand  = [ demands[i,"vu"] for i in sections]


#plot center around x = 0 
@show n = length(tendon_profile)
xmax = div(n,2)*500
@show res = mod(n+1,2)*250
@show x_range = -xmax+res:500:xmax-res



f1 = Figure(size = (1200,600))
g = f1[1,1] = GridLayout()

axs_design = Axis(g[1,1],title = "Element $element_number", titlesize = 20,
aspect = DataAspect(),
limits = (-xmax-100, xmax+100, -1.3*L, 0.6*L),
yticks = L:-100:-L,
# yminorticks = IntervalsBetween(2),
# yminorgridvisible = true,
ylabel = "y"
)

poly!(Rect( -xmax+res, -L, (n-1)*500, L*1.5), color = (:grey,0.2))
tendon = lines!(axs_design, x_range, -tendon_profile)

axs_axial  = Axis(g[2,1],aspect = 10,
limits = (-xmax, xmax, nothing, nothing),ylabel = "Axial [kN]", 
)

axs_moment = Axis(g[3,1],aspect = 10,
limits = (-xmax, xmax, nothing, nothing),ylabel = "Moment [kNm]",
)
axs_shear  = Axis(g[4,1],aspect = 10,
limits = (-xmax, xmax, nothing, nothing),ylabel = "Shear [kN]", 
)
hidexdecorations!(axs_design, grid = false)
hidexdecorations!(axs_axial, grid = false)
hidexdecorations!(axs_moment, grid = false)

lines!(axs_axial ,x_range, axial_capacity , color = :red)
lines!(axs_moment,x_range, moment_capacity, color = :blue)
lines!(axs_shear ,x_range, shear_capacity , color = :green)

lines!(axs_axial ,x_range, axial_demand ,linestyle = :dash, color = :red)
lines!(axs_moment,x_range, moment_demand,linestyle = :dash, color = :blue)
lines!(axs_shear ,x_range, shear_demand ,linestyle = :dash, color = :green)


# for (i, label) in enumerate(["Axial [kN]", "Moment [kNm]", "Shear [kN]"])
#     Box(g[i, 2], color = :gray90)
#     Label(g[i,2], label, rotation = pi/2, tellheight = false)
# end

rowgap!(g, 10)

@show yspace = maximum(tight_yticklabel_spacing!, [axs_axial, axs_shear, axs_moment])+10

axs_axial.yticklabelspace = yspace
axs_moment.yticklabelspace = yspace
axs_shear.yticklabelspace = yspace

f1

return f1
end


# ╔═╡ 5e0c7fc4-e0f8-4fb5-af21-9d71c58abdfe
designs["11"]

# ╔═╡ f7372b69-949d-437a-a38e-b4ebec8fd18c
plot_element(11, designs)

# ╔═╡ 5affd323-7913-41f2-b289-06d9c8a00d29
for i in 1:new_max_e_idx
    f = plot_element(i, designs)
    save("Results3/$i.png",f)
end

# ╔═╡ 816b861b-4768-49f8-87e4-487ea7fcd600
begin
    using CSV
    using DataFrames
    using JSON
    using Dates
    using ProgressBars
    using UnPack
    using Makie, GLMakie

end

# ╔═╡ e1da19e7-c627-4306-a0ae-f9eae259e762
begin
    println(demands[1:14,:])
    all_demands = demands
    demands = demands[1:14,:]
end

# ╔═╡ fee0e395-f598-4f24-b66b-fe25b9f3b354
begin
    #load demands into a dictionary
    demand_path = joinpath(@__DIR__, "Demands/test_input_CISBAT_dataset.json");
    open(demand_path, "r") do f
        global demands = DataFrame(JSON.parse(f, dicttype=Dict{String,Any}))
        ns = size(demands)[1]
        demands[!,:idx] = 1:ns
        println("Demands were loads from:\n", demand_path)
    end
    println(demands[1:10,:])
end

# ╔═╡ fff61fcc-5e61-4e0e-8978-689c25280bfc
begin
    using Makie, GLMakie, CairoMakie
    using JSON
    using DataFrames, CSV
end

# ╔═╡ Cell order:
# ╟─8e24d686-a754-40e3-85fe-39acd8e94ae6
# ╠═0b2b7afa-f90b-4367-9798-1348128ddc5b
# ╟─0faeda67-c89e-49fe-bb4e-d5da4684b9b4
# ╠═8b6bfef7-cf3f-4ff0-8c55-f1ac29cf1525
# ╠═816b861b-4768-49f8-87e4-487ea7fcd600
# ╟─d45ec2fe-84ea-4df6-a91a-d2e823249aad
# ╟─378a9fb6-d030-40c5-8e67-8a3e8b17697b
# ╠═e2e244f6-ec5e-42e5-b3b6-e0fabad60630
# ╠═4f705574-ec7b-4243-ae03-1da701edac30
# ╟─dc282072-4a34-485f-89f9-cc250ad93e6e
# ╠═fb79a319-57e8-4eec-b8c3-e392983fed7c
# ╟─4e048990-86f5-4bb1-86d3-7f04e5a81d8b
# ╠═fee0e395-f598-4f24-b66b-fe25b9f3b354
# ╟─748c0454-6c73-4235-b957-92d55a29001e
# ╠═a15f7dd8-5711-46ad-a01a-8911cc17aa07
# ╠═b75a66bc-ae7e-446a-83e0-66cfd1d3fdd4
# ╠═6572537d-7159-4bf2-9053-789327d7dad7
# ╠═e1da19e7-c627-4306-a0ae-f9eae259e762
# ╠═754f5262-6982-4779-991b-ed592ddc4358
# ╟─3db397c5-dea4-4bf1-b431-e2af9f6ba53d
# ╠═17db0fd9-3260-46ed-a7b6-26e2a11225c4
# ╟─d5e7772a-c241-4e4a-8bc0-5d5ebefddb1f
# ╠═1f7cfa53-4a8a-4d3c-932c-6c3d31a45d65
# ╠═70f74da6-d4b2-43aa-b3cf-be8aceff4220
# ╟─a1e47f88-ebd0-4773-88df-1ca6fc729e9c
# ╠═1c6d6600-88b6-4dfc-980e-deb6e18bb0c8
# ╟─f65446a9-ea2c-4443-9aff-354c61452bd5
# ╠═ea1d4559-5f72-4c66-b2b3-62babc8cf708
# ╠═af95e84f-c3c1-41c4-b439-e9a7ee6475d8
# ╠═54e76ba0-6e02-4f51-ba4a-a21838f6c495
# ╟─725672e0-7ae6-43da-90f8-98b3cc1911f2
# ╠═b49ed097-8ef7-4187-bfd9-f7a4400c5e30
# ╟─64a9cf05-2327-41c4-bc9f-19ff33e20049
# ╟─2ec7f3d1-9d2d-4ee9-9986-1218bbd22839
# ╠═fff61fcc-5e61-4e0e-8978-689c25280bfc
# ╟─8d6daa5d-ea86-45f1-bb0c-a12d13d34618
# ╠═d548ba9d-47ec-4f81-bb0e-94d9eb1843b9
# ╠═8eed49de-fd5e-45ad-8bab-3a590bd1b7bd
# ╠═5c26ac69-fd92-4b24-899e-890b860d9ade
# ╠═a7c95112-1348-46ef-b973-72ebcfb3276a
# ╠═5e0c7fc4-e0f8-4fb5-af21-9d71c58abdfe
# ╠═f7372b69-949d-437a-a38e-b4ebec8fd18c
# ╠═5affd323-7913-41f2-b289-06d9c8a00d29
