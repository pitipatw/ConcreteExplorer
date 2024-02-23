### A Pluto.jl notebook ###
# v0.19.38

using Markdown
using InteractiveUtils

# ╔═╡ c1f919f6-5e0f-4f2a-9fbe-c7ae52b18c32
begin
	println(pwd())
	using Pkg
	# Pkg.activate("..")
	using CSV
	using DataFrames
	using JSON
	using Dates
	using ProgressBars
	using UnPack
	using Makie, GLMakie
	using AsapSections
end


# ╔═╡ e7416f5a-d161-46db-acd7-e9e69d822abf
# ╠═╡ show_logs = false
begin
	include("Functions/definition.jl");
	include("Functions/functions.jl");
	include("Functions/structuralelement.jl")
	include("Functions/generalfunctions.jl")

	include("Functions/get_Deflection.jl")
	include("Functions/interpolations.jl");
end;

# ╔═╡ 98b2e2c0-c506-11ee-3000-a1f509a4a1a3
md"""
# PixelFrame Design Demo
"""

# ╔═╡ 650cf5aa-1db9-4f28-8325-99ae27980e58
md"""
###### To do list
1. Embodied carbon on the actual tendon profile
2. Turn everything into a gradient based optimization
"""

# ╔═╡ 02eb9acf-f1f2-496d-9784-b040ec96a4b3
md"""
Set up cell widths
"""

# ╔═╡ cfa1c5c7-3e01-40eb-9f56-e2169314c759
html"""
<style>main {max-width: 1200px;}</style>
<style>pluto-output.scroll_y {max-height: 2000px;}</style>
"""

# ╔═╡ 716f779a-9ba3-44fa-a024-89e32eb8be0a
md"""
Load required packages
"""

# ╔═╡ 41fa7984-11cb-4da2-8772-bb1e4ca40a54
md"""
Load the design catalog
Currently using version FEB6_4"""

begin
	catalog = CSV.read("src/Catalogs/FEB23_2_catalog_static.csv", DataFrame); 
	# sort!(catalog, [:carbon, :fc′, :as, :dps])
	println("The catalog was sorted by ascending order from:\ncarbon -> fc′ -> as -> dps")
	println(catalog[1:20,:])
end


# ╔═╡ 29c3fa2f-b48c-46c1-a212-48111a6db980
md"""
Load the demand points
"""

# ╔═╡ 22bd3994-0398-40c8-a2e1-594ff367810e
let
	#load demands into a dictionary
	# demand_path = joinpath(@__DIR__, "Demands/test_input_CISBAT_dataset.json");
	demand_path = joinpath(@__DIR__, "Demands/0205Full_scale_test_demands.json");
	open(demand_path, "r") do f
		global demands = DataFrame(JSON.parse(f, dicttype=Dict{String,Any}))
		ns = size(demands)[1]
		demands[!,:idx] = 1:ns
		println("There are $ns points")
		println("Demands were loads from:\n", demand_path)
	end
	println(demands[1:10,:])
end

# demands_all = demands
# 
# demands = demands[1:16, :]

# ╔═╡ d413107f-6aeb-4a94-a875-3907d240a633
md"""
From the plot, we can see that the element indices are not unique. (Different types of elements can have the same element index)
"""

# ╔═╡ f2f3ad13-fc3e-4705-8607-6c32e5d0e731
let
	println("Before Modifying the indices")
@show old_min_e_idx = minimum(demands[!, "e_idx"]);
@show old_min_s_idx = minimum(demands[!, "s_idx"]);

@show old_max_e_idx = maximum(demands[!, "e_idx"]);
@show old_max_s_idx = maximum(demands[!, "s_idx"]);


if old_min_e_idx == 0 
	demands[!,"e_idx"] .+= 1 ;
else
	println("***Element indices are not modified")
end

if old_min_s_idx == 0
	demands[!,"s_idx"] .+= 1 ;
else 
	println("***Section indices are not modified")
end

println("After Modifying the Indices")
@show new_min_e_idx = minimum(demands[!, "e_idx"]);
@show new_min_s_idx = minimum(demands[!, "s_idx"]);

@show new_max_e_idx = maximum(demands[!, "e_idx"]);
@show new_max_s_idx = maximum(demands[!, "s_idx"]);
types = unique(demands[!, :type])
global color_mapping = Dict(types .=> 1:1:length(types))
	
f_indices_check = Figure(size = (800,500))
ax_indices_check = Axis(f_indices_check[1,1], xlabel = "Global Index", ylabel = "Element Index", xticks = 0:20:250, yticks = 1:1:30)
scatter!(ax_indices_check, demands[!, "idx"], demands[!, "e_idx"], color = [color_mapping[t] for t in demands[!,:type]])
f_indices_check 
end
	

# ╔═╡ 9fb65c0f-ec4b-45ca-863c-acfe79f19c5a
md"""
Make the element indices unique
"""

# ╔═╡ 4cad732a-0ce9-4031-8cf1-d2df4c1b8502
let
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


# ╔═╡ 712343fd-57c2-4856-9f77-77ceca3c30bb
let
f_indices_check_mod = Figure(size = (1600,500))
ax_indices_check_mod = Axis(f_indices_check_mod[1,1], xlabel = "Global Index", ylabel = "Element Index", xticks = 0:20:250, yticks = 1:1:30,
title = "Modified element indices")

ax_section_indices_check_mod = Axis(f_indices_check_mod[1,2], xlabel = "Global Index", ylabel = "Element Index", xticks = 0:20:250, yticks = 1:1:30,
title = "Section indices")
ax_dps = Axis(f_indices_check_mod[1,3], xlabel = "Global Index", ylabel = "dps", xticks = 0:20:250,
title = "Section indices")
scatter!(ax_indices_check_mod, demands[!, "idx"], demands[!, "e_idx"],color = [color_mapping[t] for t in demands[!,:type]])
scatter!(ax_section_indices_check_mod, demands[!, "idx"], demands[!, "s_idx"],color = [color_mapping[t] for t in demands[!,:type]])
scatter!(ax_dps, demands[!, "idx"], demands[!, "ec_max"].*1000,color = [color_mapping[t] for t in demands[!,:type]])
f_indices_check_mod 
end
# ╔═╡ 66fea01a-ce03-4996-bd91-b6d9e91c5305
md"""
Visualize the demand points
"""

# ╔═╡ e2f306fe-7559-4856-a8b3-3757cf7fa29a
let
f1 = Figure(size = (1000,1000))
mu_max = 1.2*maximum(demands[!, :mu])
vu_max = 1.2*maximum(demands[!, :vu])

ax1 = Axis(f1[1,1], xlabel = "Moment [kNm]", ylabel = "Shear [kN]", limits = (0,mu_max,0,300), title = "Demand vs Catalog Space",
aspect = DataAspect()) 
ax2 = Axis(f1[1,2], xlabel = "Moment [kNm]", ylabel = "Shear [kN]", title= "Catalog ONLY space")
	
ax3 = Axis(f1[2,1], xlabel = "Shear [kN]" , limits = (0, vu_max, nothing, nothing), title = "Demands")
ax4 = Axis(f1[2,2], xlabel = "Shear [kN]" , limits = (0, vu_max, nothing, nothing), title = "Catalog")

ax5 = Axis(f1[3,1], xlabel = "Moment [kNm]")
ax6 = Axis(f1[3,2], xlabel = "Moment [kNm]")


scatter!(ax1, demands[!, :mu], demands[!,:vu], color = [color_mapping[t] for t in demands[!,:type]], label = demands[!, :type])
scatter!(ax1, catalog[!, :Mu], catalog[!, :Vu], colormap = :thermal , color = collect(catalog[!, :fc′]), marker = '.', alpha = 0.1, transparency = true, markersize = 5)

scatter!(ax2, catalog[!, :Mu], catalog[!, :Vu], colormap = :thermal , color = collect(catalog[!, :fc′]), marker = '.', alpha = 0.5, transparency = true, markersize = 5)



hist!(ax3, demands[!, :vu], color = :red, bins = 20)
hist!(ax4, catalog[!, :Vu], color = :blue, bins = 20)

hist!(ax5, demands[!, :mu], color = :red, bins = 20)
hist!(ax6, catalog[!, :Mu], color = :blue, bins = 20)
	
f1
end


# ╔═╡ bbbd5490-8c87-452e-976b-b6f44f91e438
md"""
Define 2 functions
1. filter_demands(demands::DataFrame, catalog::DataFrame)::Dict{Int64, Vector{Int64}}
Output => all-feasible-sections \
For each demand point, filter the feasible catalog entries and output as a dictionary in Section index => Available designs
2. find_optimum(all-feasible-sections::Dict{int64, Vector{Int64}}, demands::DataFrame})
Output => element_designs, elements_to_sections, sections_to_designs. which are Dictionaries of "first" to "second" indices.
"""

# ╔═╡ ae3344d1-89b2-42a0-b572-50bc7fc8dbd9
"""
	filter_demands!(demands::DataFrame, catalog::DataFrame)::Dict{Int64, Vector{Int64}}
For each demand point, get the feasible configuration from the catalog.
"""
function filter_demands!(demands::DataFrame, catalog::DataFrame)::Dict{Int64, Vector{Int64}}

    ns = size(demands)[1]           #total number of sections
    ne = unique(demands[!, :e_idx]) #total number of elements
    nc = size(catalog, 1)           #total number of available choices

    all_feasible_sections = Dict{Int64,Vector{Int64}}() #map between demands and indices of the feasible section.
    demands[!, "total_results"] = zeros(Int64, size(demands)[1])
    #go through each section and filter the feasible designs for the section from the catalog.
    for i = 1:ns
        en = demands[i, "e_idx"] #current element index
        sn = demands[i, "s_idx"] #current section index in that element.

        pu = demands[i, "pu"]
        mu = demands[i, "mu"]
        vu = demands[i, "vu"]
        ec_max = demands[i, "ec_max"]

        global feasible_sections = filter([:Pu, :Mu, :Vu, :dps] => (x1, x2, x3, x4) ->
                x1 >= pu &&
                x2 >= mu &&
                x3 >= vu &&
                x4 <= ec_max * 1000,
            catalog
        )
		# @show minimum(feasible_sections[:, :Mu])
        if size(feasible_sections)[1] == 0 #if the number of feasible results = 0
            println(feasible_sections[!, :ID])
            println("section $sn: element $en")
            all_feasible_sections[i] = [0]
            demands[i, "total_results"] = 0
            # println(outr)
        else
			# @show minimum(feasible_sections[!, :Mu])
            all_feasible_sections[i] = feasible_sections[!, :ID]
            demands[i, "total_results"] = length(feasible_sections[!, :ID])
            # println(outr)
        end
    end

    return all_feasible_sections
end


# ╔═╡ d02089f3-5123-4cdb-8f5f-810a393e5e4e
# begin 
# 	ne = unique(demands[!, :e_idx])
# 	elements_to_sections = Dict(k => Int[] for k in ne)
# 	for i in 1:size(demands)[1]
# 		en = demands[i, :e_idx]
# 		# sn = demands[i, "s_idx"]
# 		push!(elements_to_sections[en], demands[i, :idx])
# 	end
# 	println("Elements -> Sections in them")
# 	for i in 1:size(unique(demands[!, :e_idx]))[1]
# 		println(i, " => ",elements_to_sections[i])
# 	end
# end



# ╔═╡ 0d5119c9-8874-4efd-8400-22c359b804cf
"""
	function find_optimum(all_feasible_sections::Dict{Int64, Vector{Int64}}, demands::DataFrame)
Find the optimum result for each element.
    
    !! not optimum yet.
Constraints for the same element.


"""
function find_optimum(all_feasible_sections::Dict{Int64, Vector{Int64}}, demands::DataFrame)
	total_number_of_sections = size(demands)[1] #get total number of section points.
    ne = unique(demands[!, :e_idx]) #list of element index
	
	#Outputs
	elements_to_sections = Dict(k => Int[] for k in ne) #Map element indices to indices of sections on those elements.
    elements_designs = Dict(k => [[]] for k in unique(demands[!, :e_idx]))  #element index to list of designs
    sections_to_designs = Dict(k => Vector{Float64}() for k in demands[!, :idx]) #each section to its design

	for i in 1:total_number_of_sections
		e = demands[i, :e_idx] 
		push!(elements_to_sections[e], demands[i, :idx])
	end
	
	println("Elements -> Sections in them")
	for i in 1:size(ne)[1]
		println(i, " => ",elements_to_sections[i])
	end

    #find as, and fpe that appear in all sections in an element.
    for i in ne #go through each element.
        println("Element $i out of $(length(ne)) elements")
		
        sections = elements_to_sections[i] #sections associated with this element
        ns = length(sections)     # number of sections in this element 

        #start from the middle-ish section (n/2 or (n-1)/2), it must be the section with the smallest available configuration.
        #note that section is in the form of 1,2,3,..., ns.
        mid = div(ns, 2)
        
        #get the feasible designs for the middle section
        feasible_idx = all_feasible_sections[sections[mid]]
		println(size(feasible_idx)[1], " available sections")

        #catalog was already sorted, so I think we can leave this part, just filter, to save time.
        mid_catalog = sort(catalog[feasible_idx, :], [:carbon, :fc′, :dps])

        #now, loop each design in the sub catalog, see if "as" and "fpe" are available in all sections.
        #if not, remove that design from the sub catalog.
        #if yes, keep it.

        #select each design, check if as and fpe exist for the the section
        global final_d_idx = 0 
		total_mid_catalog = size(mid_catalog)[1]
		global found_all = true #make it a global variable.
        for d_idx in 1:total_mid_catalog # go through every possible mid catalog.
            d = mid_catalog[d_idx, :]
			found_all = true
            # all_as  = true
            # all_fpe = true
            serviceability_check = true
            for s in sections #check if as and fpe occurs in other feasible designs of other sections.
				println("Check section $s")
                #if not, go to the next choice.

                #If matching fc′, make sure that the fR1 and fR3 are also matched.

                # if !(d[:fc′] ∈ catalog[all_feasible_sections[s], :fc′])
				#also check the fR1 and fR3
                #     all_fc′ = false
                #     break
                # end
				check_as =  d[:as] ∈ catalog[all_feasible_sections[s], :as]
				check_fpe = d[:fpe] ∈ catalog[all_feasible_sections[s], :fpe]
				
				if !check_as || !check_fpe #not found, move to the next design of the mid section.
                    found_all = false
					println("Section $s fails, restarting...")
                    break
                end
            end

			if found_all
				println("Found all at element $i")
				final_d_idx = d_idx
				break
			end
        end

		#If the loop finishes without false, that's a hit! Otherwise, not found
		if !found_all
			println("Warning, can't find the solution for element $i")
			#output empty parameters (In this case 0)
			elements_designs[i] = [0]
    		sections_to_designs[elements_to_sections] .= [0.]
		else
			println("making element $i")
	        #get the first one, they will appear in the entire thing anyway.
	        # this_fc′ = mid_catalog[global_d, :fc′]
			# this_fR1 = mid_catalog[global_d, :fR1]
			# this_fR3 = mid_catalog[global_d, :fR3]
			println("final_d_idx is $final_d_idx")
	        this_fpe = mid_catalog[final_d_idx, :fpe] 
	        this_as = mid_catalog[final_d_idx, :as]
	
	        sections_designs = Vector{Vector}(undef, ns)
	        for is in eachindex(elements_to_sections[i])
	            #current section index
	            s = elements_to_sections[i][is]
	
	            feasible_idx = all_feasible_sections[s] # all feasible sections for this section.
				# println("Feasible Catalog")
	            # fc′_fpe_as(fc′::Float64, fpe::Float64, as::Float64) = fc′ == this_fc′ && fpe == this_fpe && as == this_as
	            fpe_as(fpe::Float64, as::Float64) = fpe == this_fpe && as == this_as
	
	            # this_catalog = filter([:fc′, :fpe, :as] => fc′_fpe_as, catalog[output_results[s], :])
	            this_catalog = filter([:fpe, :as] => fpe_as, catalog[feasible_idx, :])
	
	            sort!(this_catalog, [:carbon, :dps]) #the lowest carbon will be the first index.
	            select_ID = this_catalog[1, :ID] #The first one is the lowest.
	            sections_designs[is] = collect(catalog[select_ID, :])
	        end

			# #Create a PixelFrame element -> Find the deflection of this element. (Beam, Column, etc).
	  #       #Create a pixelframeelement and/or section here with the given parameters 
	  #       L, t, Lc = [205.0 35.0 30.0] #Should make this tie to the catalog, but for now we only have 1 configuration of L,t,Lc.
	  #       compoundsection =  make_Y_layup_section(L, t, Lc)
	  #       pixelframeelement = PixelFrameElement() 
			# #Modified the fields inside pixelframeelement so they are met with the current configuration.
			# pixelframeelement.compoundsection = compoundsection
	  #       Le = ns*500.0 #500 mm per section.
	
	  #       #could do a case where input only the variables -> the dependent variables come later.
	  #       pixelframeelement.fc′ = this_fc′ # Concrete strength [MPa] ****Should update on the test day using cylinder test***
	  #       # Ec = 4700.0*sqrt(fc′) # MPa  ACI fc-> Concrete modulus relationship [MPa]
	  #       pixelframeelement.Ec = 58000.0 # MPa  from the cylinder test
	  #       pixelframeelement.Eps = 70000.0 #Post tensioning steel modulus [MPa]
	  #       pixelframeelement.fpy = 0.002 * pixelframeelement.Eps #MPa  
	  #       #Safe load on the website https://www.engineeringtoolbox.com/wire-rope-strength-d_1518.html 
	  #       # is ~ 150 MPa. Currently 140 MPa :)
	
	  #       # PixelFrame section/element properties
	  #       centroid_to_top = 100.0 #[mm] ~half of 205mm
	  #       pixelframeelement.em = mid_catalog[global_d, :dps]# Eccentricity at the middle of the member [mm]
	  #       pixelframeelement.es = 0.0 # Eccentricity at the support of the member   [mm]
	  #       pixelframeelement.em0 = mid_catalog[global_d, :dps] # Initial eccentricity at the midspan        [mm]
			# # Initial distance from the top to the point of application of the 		load [mm]
	  #       pixelframeelement.dps0 = centroid_to_top + pixelframeelement.em0 
	  #       pixelframeelement.Ls = Le/4 # Distance from support to the first load point [mm]
	  #       pixelframeelement.Ld = Le/4 # Distance from support to the first deviator [mm]
	  #       pixelframeelement.L = Le # Total length of the member [mm]
	  #       # two 1/4" bars with 1200 lb capacity
	  #       pixelframeelement.Aps = this_as # Total area of the post tensioned steel [mm2]
	  #       #Pure concrete area = 18537.69 mm2
	  #       #Transformed steel area = 347.96 mm2 
			# # Transformed area of the cross section [mm2] (= Concrete area if there is no embedded rebars)
	  #       pixelframeelement.Atr = compoundsection.area 
	  #       pixelframeelement.Itr = compoundsection.Ix #moment of inertia [mm4], no embedded steel, therefore, only from concrete.
	  #       # pixelframeelement.Itr = 1.082e+8 #this number includes deviated steels.
		 #    # Elastic modulus of the concrete section from the centroid to extreme tension fiber [mm3]
	  #       pixelframeelement.Zb = pixelframeelement.Itr/centroid_to_top 
	  #       # If there are multiple materials, transformed section geometry is needed for Zb (and everything related to section area)
	
	  #       #forces
	  #       pixelframeelement.w = pixelframeelement.Atr / 10^9 * 2400.0 * 9.81 # Selfweight [N/mm]
	  #       pixelframeelement.mg = pixelframeelement.w * pixelframeelement.L^2 / 8.0 # Moment due to selfweight [Nmm]
	  #       pixelframeelement.fr = 0.7 * sqrt(pixelframeelement.fc′) # Concrete cracking strenght [MPa]
	  #       pixelframeelement.r = sqrt(pixelframeelement.Itr / pixelframeelement.Atr) # Radius of gyration [mm]
	  #       pixelframeelement.ps_force = pixelframeelment.Aps*this_fpe # Post tensioning force [N]
	  #       pixelframeelement.Mdec = pixelframeelement.ps_force*pixelframeelement.em
	  #       pixelframeelement.concrete_force = pixelframeelement.ps_force*cos(24.0*pi/180.0) # should use actual value 
	  #       pixelframeelement.fpe = pixelframeelement.ps_force/pixelframeelement.Aps # Effective post tensioning stress [MPa] 
	  #       pixelframeelement.ϵpe = pixelframeelement.fpe / pixelframeelement.Eps # Effective post tensioning strain [mm/mm]
	  #       #find moment due to the applied force.
	  #       pixelframeelement.ϵce = pixelframeelement.ps_force*pixelframeelement.em/pixelframeelement.Zb/pixelframeelement.Ec - pixelframeelement.concrete_force/pixelframeelement.Atr/pixelframeelement.Ec # effetive strain in the concrete [mm/mm]
	  #       #for using test setup
	  #       pixelframeelement.test = false #This is not the half-scale test beam anymore.
	
	  #       #for the current stage of the model,
	  #       #we will have 2 deviators, with interpolated 1/4 of the length of the element.
	  #       #with this, <Need proved, but this should be more conservative???>. Less post tensioned than the actual beam.
	  #       #Current PixelFrameSize
	        
	  #       load_m = 4*demands[mid, :mu]/Le*1000 #N
	  #       load_v = demands[mid, :vu]*1000
	
	  #       load = load_m > load_v ? load_m : load_v
	  #       # @show load 
	  #       dis_history, P =  get_Deflection(pixelframeelement, load, loadstep = 1000)
			# @show length(dis_history)
	  #       δ = dis_history[end]
	
	        # elements_designs[i] = [sections_designs, δ, δ/(Le/240)]
			elements_designs[i] = sections_designs #, δ, δ/(Le/240)]

	
	    end
		
	end
	#at this point, we have all of the designs of the elements!!!
		# println(typeof(sections_to_designs))
	    for i in ne #loop each element
			#i is an index of an element.
	        @show sections = elements_to_sections[i] #sections numbers in that element.
			#maximum sections 
			# max_section = maximum(sections)
			@show elements_designs[i]
	        for design_idx in eachindex(sections)
				# @show typeof(elements_designs[i])
				# @show typeof(elements_designs[i][design_idx])
				# println(elements_designs[i][design_idx])

				# @show sections[design_idx]
				# @show sections_to_designs[sections[design_idx]]
	            sections_to_designs[sections[design_idx]] = elements_designs[i][design_idx]
	        end
	    end
    return elements_designs, elements_to_sections, sections_to_designs
end

# ╔═╡ df8a7202-49c9-4e8a-81eb-aa63cc410888
let
global all_feasible_sections= filter_demands!(demands,catalog)
end

#select a section to see the available designs
section_number = 1
figure_check_section = Figure(size = (500,500))
ax_section = Axis(figure_check_section[1,1], xlabel = "Moment [kNm]", ylabel = "Shear [kN]", title = string(section_number))
for d in 1:length(all_feasible_sections[section_number])
	scatter!(ax_section, catalog[d,:Mu], catalog[d,:Vu])
end

for i in 1:length(all_feasible_sections)
	println("Section $i")
	println(length(all_feasible_sections[i]))
end

# ╔═╡ 0176e9c5-de48-4163-9cfb-b909177869dc
# println(all_feasible_sections)

# ╔═╡ 5df9ac50-bdef-4cfb-845f-8d84b9be1ef8
elements_designs, elements_to_sections, sections_to_designs  = find_optimum(all_feasible_sections, demands)
println(elements_designs)


elements_designs_fielded = Vector{Dict{String,Real}}()
# for i in eachindex(elements_designs)
open("src/Results/designs_results_23_02.json","w") do f
    JSON.print(f, elements_designs)
end


open("src/Results/sections_to_designs_23_02.json","w") do f
    JSON.print(f, sections_to_designs)
end

using Makie, GLMakie, CairoMakie
using JSON
using DataFrames, CSV

designs = JSON.parsefile(joinpath(@__DIR__,"Results/designs_results_06_02.json"), dicttype = Dict{String,Vector{Vector{Float64}}});
mapping_strings = ["ID","fc′", "fR1", "fR3", "as" ,"dps", "fpe", "Pu" ,"Mu", "Vu,", "carbon", "L", "t", "Lc","T", "catalog_id"]
for i in 1:length(elements_designs)
	for s in eachindex(elements_designs[i])
		global_s_index = elements_to_sections[i][s]
		println(global_s_index)
		# @show  Dict(mapping_strings .=> elements_designs[parse(Int64,i)][s])
		# @show elements_designs_fielded[global_s_index]
		push!(elements_designs_fielded, Dict(mapping_strings .=> vcat(global_s_index-1, elements_designs[i][s])))
		# elements_designs_fielded[global_s_index] =
	end
end 

open("src/Results/09_02_designs_results_fielded.json","w") do f
    JSON.print(f, elements_designs_fielded)
end

ne = length(designs)
println("There are $ne elements.")

function plot_element(element_number::Int64, designs::Dict;L::Float64 = 250.0)


	sections = elements_to_sections[element_number]
	element_number = string(element_number)
	
	L = 205
	tendon_profile = [i[5] for i in designs[element_number]]
	axial_capacity  = [i[7] for i in designs[element_number]]
	moment_capacity = [i[8] for i in designs[element_number]]
	shear_capacity = [i[9] for i in designs[element_number]]
	
	
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
	

	#test 
	plot_element(1, designs)
	
	for i in 1:19
		f = plot_element(i, designs)
		save("src/Results5/$i.png",f)
	end

	#summarize the result. 
	function get_design_properties(sections_to_designs::Dict{Int64, Vector{Float64}}, idx::Int64)
		output = Vector{Float64}(undef, length(sections_to_designs))
		for i in eachindex(sections_to_designs)
			@show i
			output[i] = sections_to_designs[i][idx]
		end
		return output
	end


	
	all_fc′  = get_design_properties(sections_to_designs,1)
	all_fR1  = get_design_properties(sections_to_designs,2)
	all_fR3  = get_design_properties(sections_to_designs,3)
	all_as   = get_design_properties(sections_to_designs,4)
	all_dps  = get_design_properties(sections_to_designs,5)
	all_fpe  = get_design_properties(sections_to_designs,6)
	stack_name = hcat(string.(all_fc′,"_", all_fR1,"_", all_fR3))
	@show unique_stack_name = unique(stack_name)
	MK_file_prep = DataFrame(:fc′ => all_fc′, :fR1 => all_fR1, :fR3=>all_fR3)

	csv_fc′ = Vector{Float64}()
	csv_fR1 = Vector{Float64}()
	csv_fR3 = Vector{Float64}()

	for i in 1:length(unique_stack_name)
		@show vals = parse.(Float64,split(unique_stack_name[i], "_"))
		@show typeof(vals[1])
		@show typeof(vals[2])
		@show typeof(vals[3])
		push!(csv_fc′, vals[1])
		push!(csv_fR1, vals[2])
		push!(csv_fR3, vals[3])
	end
	
	csv_output = DataFrame(:fc′=> csv_fc′, :fR1=> csv_fR1, :fR3=>csv_fR3)
	CSV.write(joinpath(@__DIR__, "mix_specs.csv"), csv_output)
	
	
	
	#each pair, plots them dots and x and a line connecting them together. 




	f_final = Figure(size= (500,500))
	ax1 = Axis(f_final[1,1], xlabel = "Moment [kNm]", ylabel = "Shear [kN]", title = "Demands vs Designs")
	demand_points = hcat(demands[!, :mu] , demands[!,:vu])
	design_points = hcat(get_design_properties(sections_to_designs,8),get_design_properties(sections_to_designs,9))
	for i in 1:size(demand_points)[1]
		x1 = demand_points[i,1]
		y1 = demand_points[i,2]
		x2 = design_points[i,1]
		y2 = design_points[i,2]
		@assert x2>x1
		@assert y2>y1
		u = x2-x1
		v = y2-y1
		arrows!([x1],[y1],[u],[v], arrowsize = 5)
	end
	scatter!(ax1, demand_points[:,1],demand_points[:,2], color = :red, markersize = 10)
	scatter!(ax1, design_points[:,1],design_points[:,2], color = all_fc′, market_size = 10)

	f_final



	f = Figure(size = (800, 800))
	Axis(f[1, 1], backgroundcolor = "black")
	
	xs = LinRange(0, 2pi, 20)
	ys = LinRange(0, 3pi, 20)
	us = [sin(x) * cos(y) for x in xs, y in ys]
	vs = [-cos(x) * sin(y) for x in xs, y in ys]
	strength = vec(sqrt.(us .^ 2 .+ vs .^ 2))
	
	arrows!(xs, ys, us, vs, arrowsize = 10, lengthscale = 0.3,
		arrowcolor = strength, linecolor = strength)
	
	f

# ╔═╡ Cell order:
# ╟─98b2e2c0-c506-11ee-3000-a1f509a4a1a3
# ╟─650cf5aa-1db9-4f28-8325-99ae27980e58
# ╟─02eb9acf-f1f2-496d-9784-b040ec96a4b3
# ╠═cfa1c5c7-3e01-40eb-9f56-e2169314c759
# ╟─716f779a-9ba3-44fa-a024-89e32eb8be0a
# ╠═c1f919f6-5e0f-4f2a-9fbe-c7ae52b18c32
# ╠═e7416f5a-d161-46db-acd7-e9e69d822abf
# ╟─41fa7984-11cb-4da2-8772-bb1e4ca40a54
# ╠═f074d05a-4a52-485f-b0cb-4a79a6fa42b3
# ╟─29c3fa2f-b48c-46c1-a212-48111a6db980
# ╠═22bd3994-0398-40c8-a2e1-594ff367810e
# ╟─d413107f-6aeb-4a94-a875-3907d240a633
# ╠═f2f3ad13-fc3e-4705-8607-6c32e5d0e731
# ╟─9fb65c0f-ec4b-45ca-863c-acfe79f19c5a
# ╠═4cad732a-0ce9-4031-8cf1-d2df4c1b8502
# ╠═712343fd-57c2-4856-9f77-77ceca3c30bb
# ╟─66fea01a-ce03-4996-bd91-b6d9e91c5305
# ╠═e2f306fe-7559-4856-a8b3-3757cf7fa29a
# ╟─bbbd5490-8c87-452e-976b-b6f44f91e438
# ╠═ae3344d1-89b2-42a0-b572-50bc7fc8dbd9
# ╠═d02089f3-5123-4cdb-8f5f-810a393e5e4e
# ╠═0d5119c9-8874-4efd-8400-22c359b804cf
# ╠═df8a7202-49c9-4e8a-81eb-aa63cc410888
# ╠═0176e9c5-de48-4163-9cfb-b909177869dc
# ╠═5df9ac50-bdef-4cfb-845f-8d84b9be1ef8
