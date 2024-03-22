"""
Define 2 functions
1. filter_demands(demands::DataFrame, catalog::DataFrame)::Dict{Int64, Vector{Int64}}
Output => all-feasible-sections \
For each demand point, filter the feasible catalog entries and output as a dictionary in Section index => Available designs
2. find_optimum(all-feasible-sections::Dict{int64, Vector{Int64}}, demands::DataFrame})
Output => element_designs, elements_to_sections, sections_to_designs. which are Dictionaries of "first" to "second" indices.
"""

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

        pu = abs(demands[i, "pu"])
        mu = demands[i, "mu"]
        vu = demands[i, "vu"]
        ec_max = demands[i, "ec_max"]*1000
		T_string = demands[i, "type"]
		if T_string == "primary" || T_string == "secondary"
			T = 3	
		elseif T_string == "columns"
			T = 2
		else 
			println("Warning, invalid type")
		end

        feasible_sections = filter([:Pu, :Mu, :Vu, :dps, :T] => (x1, x2, x3, x4, x5) ->
                x1 >= pu &&
                x2 >= mu &&
                x3 >= vu &&
                x4 <= ec_max &&
				mod(x5,T) == 0, # T is the numbers, 3 belongs ot Primary and Secondary, 2 belongs to Columns, which has 2 and 4. 
            catalog
        )
        @assert minimum(feasible_sections[!, :Vu]) >= vu
		@assert minimum(feasible_sections[:, :Mu]) >= mu
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


"""
	function find_optimum(all_feasible_sections::Dict{Int64, Vector{Int64}}, demands::DataFrame)
Find the optimum result for each element.
Design map from element number -> sections -> designs
    !! not optimum yet.
Constraints for the same element.
"""
function search_design(all_feasible_sections::Dict{Int64, Vector{Int64}}, demands::DataFrame;
	demands_ser::DataFrame = DataFrame(),
	start_mid::Bool = true)

	if size(demands_ser)[1] == 1 
		println("Couldn't find service demands...")
		println("Using ultimate demands as service demands (conservative)")
		demands_ser = demands 
	end

	total_number_of_sections = size(demands)[1] #get total number of section points.
    ne = unique(demands[!, :e_idx]) #list of element index
	
	# Pre allocate outputs
	elements_to_sections = Dict(k => Int[] for k in ne) #Map element indices to indices of sections on those elements.
    elements_designs = Dict(k => Dict() for k in unique(demands[!, :e_idx]))  #element index to list of designs
    sections_to_designs = Dict(k => Vector{Float64}() for k in demands[!, :idx]) #each section to its design
	sections_to_designs = Dict(k => Dict() for k in demands[!, :idx]) 
	skipped_elements = Vector{Int64}()

	catalog_keys::Vector{Symbol} = [:fc′, :dosage, :fR1, :fR3, :as, :dps, :fpe, :fps, :Pu, :Mu, :Vu, :carbon, :L, :t, :Lc, :T, :ID];
	
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
		this_type = demands[sections[1],:type]

		#check if there is a section in that element that's not feasible (available section -> [0]),
		# if this exists, put result into a vector of 0 and skip this element. 
		for s in sections
			if all_feasible_sections[s][1] == 0 
				push!(skipped_elements, i)
				break
			end
		end

		if i ∈ skipped_elements
			elements_designs[i] = Dict(sections .=> 0)
			println("Element $i skipped")
			continue
		end

		if start_mid
        #start from the middle-ish section (n/2 or (n-1)/2), it must be the section with the smallest available configuration for moment.
        #note that section is in the form of 1,2,3,..., ns.
        	mid = div(ns, 2)
		else 
			#otherwise, starts from 1 (shear controlled?)
			mid = 1
		end
        
        #get the feasible designs for the middle section
        feasible_idx = all_feasible_sections[sections[mid]]
		println(size(feasible_idx)[1], " available sections")

        #catalog was already sorted, so I think we can leave this part, just filter, to save time.
        mid_catalog = sort(catalog[feasible_idx, :], [:carbon, :fc′, :fpe, order(:dps, rev = true)])

        #now, loop each design in the sub catalog, see if "as" and "fpe" are available in all sections.
        #if not, remove that design from the sub catalog.
        #if yes, keep it.

        #select each design, check if as and fpe exist for the the section
		total_mid_catalog = size(mid_catalog)[1]
		final_design_index = 0 
		found_all = true #define here so it is availables outside of the for loop.
        for d_idx in 1:total_mid_catalog # go through every possible mid catalog.

            current_mid_design = mid_catalog[d_idx, :]
			this_fpe  = current_mid_design[:fpe]
			this_as   = current_mid_design[:as]
			this_type = current_mid_design[:T]
			this_dps = current_mid_design[:dps]

			this_L = current_mid_design[:L]
			this_t = current_mid_design[:t]
			this_Lc = current_mid_design[:Lc]

			#create a filter function that check every constrained parameter at once.
			found_all = true

			if this_type == "secondary"
				#if it's a secondary beam, additionally fix the dps.
				fpe_as_dps_type(fpe::Float64, as::Float64, dps::Float64,type::Float64, L::Real, t::Real, Lc::Real) = fpe == this_fpe && as == this_as && dps == this_dps && type == this_type && L == this_L && t == this_t && Lc == this_Lc
			else
				fpe_as_type(fpe::Float64, as::Float64, type::Float64, L::Real, t::Real, Lc::Real) = fpe == this_fpe && as == this_as && type == this_type && L == this_L && t == this_t && Lc == this_Lc
			end
            # serviceability_check = true # Will add this.

            for s in sections #check if as and fpe occurs in other feasible designs of other sections.
				println("Check section $s")
				feasible_sections = all_feasible_sections[s]
				section_feasible_catalog = catalog[feasible_sections,:]
      
				#a check function, ensures that these parameters happen at the same time.
				if this_type == "secondary" 
					this_catalog = filter([:fpe, :as, :dps,:T, :L,:t,:Lc] => fpe_as_dps_type, section_feasible_catalog)
				else
					this_catalog = filter([:fpe, :as, :T, :L,:t,:Lc] => fpe_as_type, section_feasible_catalog)
				end

				if size(this_catalog)[1] == 0
                    found_all = false #if this happens at the last possible mid section, the found_all would trigger the next if-else check on the next block.
					println("Section $s fails, restarting...")
                    break
                end
            end

			if found_all
				println("Found all at element $i")
				final_design_index = d_idx
				break
			end
        end

		#If the loop finishes without false, that's a hit! Otherwise, not found
		if !found_all
			println("!!!!!Warning, couldn't find the solution for element $i")
			println("Outputting empty parameters (In this case 0)")
			elements_designs[i] = [0]
    		sections_to_designs[elements_to_sections][catalog_keys] .= [0.]
		else #we start the searching here.
			println("###")
			println("Making element $i")
			println("final_design_index is $final_design_index")
			#This could be removed since the final_design_index should be the same as the d_index defined in the loop.
	        this_fpe  = mid_catalog[final_design_index, :fpe] 
	        this_as   = mid_catalog[final_design_index, :as]
			this_type = mid_catalog[final_design_index, :T]
			this_dps = mid_catalog[final_design_index, :dps]
			this_L = mid_catalog[final_design_index, :L]
			this_t = mid_catalog[final_design_index, :t]
			this_Lc = mid_catalog[final_design_index, :Lc]

			sections_designs = Dict{Int64,  Dict{Symbol,Union{Float64,Int64,String}}}(k =>Dict() for k in 1:ns) 
	        for is in eachindex(elements_to_sections[i])

	            #current section index
	            s = elements_to_sections[i][is]
	
	            feasible_idx = all_feasible_sections[s] # all feasible sections for this section.
				# fpe_as_type(fpe::Float64, as::Float64, type::Float64) = fpe == this_fpe && as == this_as && type == this_type
				section_feasible_catalog = catalog[feasible_idx,:]
	            # this_catalog = filter([:fc′, :fpe, :as] => fc′_fpe_as, catalog[output_results[s], :])
				if this_type == "secondary"
					fpe_as_dps_type(fpe::Float64, as::Float64, dps::Float64,type::Float64, L::Real, t::Real, Lc::Real) = fpe == this_fpe && as == this_as && dps == this_dps && type == this_type && L == this_L && t == this_t && Lc == this_Lc
					this_catalog = filter([:fpe, :as, :dps,:T, :L,:t,:Lc] => fpe_as_dps_type, section_feasible_catalog)
					sort!(this_catalog, [:carbon, order(:dps, rev=true)] ) #the lowest carbon then, dps will be the first index.

				else
					fpe_as_type(fpe::Float64, as::Float64, type::Float64, L::Real, t::Real, Lc::Real) = fpe == this_fpe && as == this_as && type == this_type && L == this_L && t == this_t && Lc == this_Lc
					this_catalog = filter([:fpe, :as, :T, :L,:t,:Lc] => fpe_as_type, section_feasible_catalog)
					# sort!(this_catalog, [:carbon, order(:dps, rev=true)] ) #should pick the lowest dps first, so it matched with the demands
					#otherwise, it would be confusing on why the capacity is so high.
					sort!(this_catalog, [:carbon, :dps]) #should pick the lowest dps first, so it matched with the demands

				end
	            # this_catalog = filter([:fpe, :as, :T] => fpe_as_type, catalog[feasible_idx, :])
				if size(this_catalog)[1] == 0 
					@show this_catalog
				end

				
	            select_ID = this_catalog[1, :ID] #The first one is the lowest embodied carbon and highest dps.

				#add maxmimum and minimum dps for this part.
				#filter again for all of the same configuration except dps.
				this_fc′ = this_catalog[1,:fc′]

				#do this so the sections in the band are still the same, but different dps.
				final_cut(fc′::Float64, fpe::Float64, as::Float64) = fc′ == this_fc′ && fpe == this_fpe && as == this_as
				constrained_catalog = filter([:fc′ ,:fpe, :as] => final_cut , catalog[feasible_idx, :])

				maximum_dps = maximum(constrained_catalog[!, :dps])
				minimum_dps = minimum(constrained_catalog[!, :dps])

				# @show vcat(collect(catalog[select_ID, :]), [maximum_dps, minimum_dps])
	            sections_designs[is] = Dict(vcat(catalog_keys, [:max_dps, :min_dps,]) .=> vcat(collect(catalog[select_ID, :]), [maximum_dps, minimum_dps]))
	        end

	  #       Create a PixelFrame element -> Find the deflection of this element. (Beam, Column, etc).
	  #		  Then, check with the limit L/240, if pass, move on to the next element, otherwise, go to the next configuration.
	  #       Create a pixelframeelement and/or section here with the given parameters 
	        # L, t, Lc = [205.0 35.0 30.0] #Should make this tie to the catalog, but for now we only have 1 configuration of L,t,Lc.
			L = mid_catalog[final_design_index, :L] 
			t = mid_catalog[final_design_index, :t] 
			Lc = mid_catalog[final_design_index, :Lc] 
			if this_type == 2.0
				compoundsection =  make_X2_layup_section(L, t, Lc)
			elseif this_type == 3.0
				compoundsection =  make_Y_layup_section(L, t, Lc)
			elseif this_type == 4.0
				compoundsection =  make_X4_layup_section(L, t, Lc)
			else
				println("Invalid type")
			end
	  #       pixelframeelement = PixelFrameElement() 
	  #		# Modified the fields inside pixelframeelement so they are met with the current configuration.
	  #		# pixelframeelement.compoundsection = compoundsection
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
	  #		# # Transformed area of the cross section [mm2] (= Concrete area if there is no embedded rebars)
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
	#   if δ> L/240 
	# continue
	#   end 
	  # elements_designs[i] = [sections_designs, δ, δ/(Le/240)]
			elements_designs[i] = sections_designs #, δ, δ/(Le/240)]
	    end
	end
	#at this point, we have all of the designs of the elements

	#Determine the axial force.
	
	# find the angle between the support face and the next deviated shape. 

	μs = 0.3
	for e in ne
		if e ∉ skipped_elements
		element_designs = elements_designs[e]
		ns = length(element_designs)
		support_dps = element_designs[1][:dps]
		next_dps = element_designs[2][:dps]
	
		another_support_dps = element_designs[ns][:dps]
		another_next_dps = element_designs[ns-1][:dps]
			
		#There is a symmetry (there must be a symmetry)
		@assert support_dps == another_support_dps "Symmetry Error at supports: $support_dps ≠ $another_support_dps."
		@assert next_dps == another_next_dps "Symmetry Error next to supports: $next_dps ≠ $another_next_dps."
		
		as = element_designs[1][:as]
		fps = element_designs[1][:fps]

		L = 500.0 # That's the distance between the 2 deviated tendons (conservative)
		θ = atan((next_dps-support_dps)/L)
		# rad2deg(θ)
		axial_component = as*fps*cos(θ)/1000 #kN
		maxV = 0
		for i in eachindex(element_designs)
			# println(element_designs[i][:Vu])
			if element_designs[i][:Vu] > maxV 
				maxV = element_designs[i][:Vu]
			end
		end

		required_normal_force = maxV/μs #kN
		additional_force = required_normal_force - axial_component
	
		for s in eachindex(element_designs)
			element_designs[s][:axial_force] = additional_force
		end
	end
	end

	#now we are generating different output format for the ease of use in the later steps.
	for e in ne #loop each element
		if e ∉ skipped_elements 
					#i is an index of an element.
		sections = elements_to_sections[e] #sections numbers in that element.

		# elements_designs[i]
		for design_idx in eachindex(sections)
			# @show elements_designs[i][design_idx]
			sections_to_designs[sections[design_idx]] = elements_designs[e][design_idx]
		end
	end
	end
    return elements_designs, elements_to_sections, sections_to_designs, skipped_elements
end