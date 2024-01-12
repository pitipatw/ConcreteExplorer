using CSV, DataFrames, JSON
using Dates

"""
catalog format
fc', as, ec, fpe, Pu, Mu, Vu, embodied
"""

catalog = CSV.read(joinpath(@__DIR__,"Catalogs/catalog_static.csv"), DataFrame);

sort!(catalog, [:carbon, :fc′, :as, :ec])
println("The catalog was sorted by ascending order from:")
println("carbon -> fc′ -> as -> ec")
println(catalog[1:100,:])

#load demands into a dictionary
demand_path = joinpath(@__DIR__, "Demands/CISBAT_test_input.json");
open(demand_path, "r") do f
    global demands = DataFrame(JSON.parse(f, dicttype=Dict{String,Any}))
    ns = size(demands)[1]
    demands[!,:idx] = 1:ns
    println("Demands were loads from:\n", demand_path)
end

#python used 0 index, here, we shift those by 1.
demands[!,"e_idx"] .+= 1
demands[!, "s_idx"] .+=1
println("Section and Element indices were shifted by 1, (0 to 1 based)")

#properly label the element index. Now there are repetitions.
#Need to generate a new version that does not repeat.
println("element indices were redundant, now they are fixed")
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

println(demands)
# sort!(demands, [:e_idx, :s_idx])

function filter_demands(demands::DataFrame, catalog::DataFrame)

#map from element idx to section idx.

ns = size(demands)[1]          #total number of sections
ne = unique(demands[!,:e_idx]) #total number of elements
nc = size(catalog,1)           #total number of available choices

output_results = Dict{Int64, Vector{Int64}}() #map between demands and indx of the feasible section they have

demands[!, "num_results"] = zeros(Int64, size(demands)[1])
#go through each section and filter the feasible designs from the catalog.
for i = 1:ns
    en = demands[i, "e_idx"]
    sn = demands[i, "s_idx"]
    # push!(elements_to_sections[en], sn)

    pu = demands[i,"pu"]
    mu = demands[i,"mu"]
    vu = demands[i,"vu"]
    ec_max = demands[i,"ec_max"]

    global feasible_sections = filter([:Pu,:Mu,:Vu,:ec] => (x1,x2,x3,x4) -> 
    x1>pu &&
    x2>mu &&
    # x3 > vu &&
    x4<= ec_max*1000, 
    catalog
    )

    if size(feasible_sections)[1] == 0 #if the number of feasible results = 0
        println(feasible_sections[!, :ID])
        println("section $sn: element $en")
        println(vu)
        output_results[i] = [1]
        demands[i, "num_results"] = 0
        # println(outr)
    else
        output_results[i] = feasible_sections[!,:ID]
        demands[i, "num_results"] = length(feasible_sections[!, :ID])
        # println(outr)
    end
end

return output_results
end

output_results = filter_demands(demands,catalog)

"""
Find the optimum result for each element.
    
    !! not optimum
For the same element, will use the same fc′ steel size and post tensioning stress.

"""
function find_optimum(output_results, demands)
ns = size(demands)[1]
ne = unique(demands[!,:e_idx]) #python starts at 0
elements_designs = Dict(k=> [[]] for k in unique(demands[!,"e_idx"])) #element index to list of designs

elements_to_sections = Dict(k=> Int[] for k in unique(demands[!,"e_idx"]))
for i = 1:ns
    en = demands[i, "e_idx"]
    # sn = demands[i, "s_idx"]
    push!(elements_to_sections[en], i)
end


#find fc′, as, and fpe that appear in all sections in an element.
for i in ne #loop each element
    println("working on element $i out of $(length(ne)) elements")
    sections = elements_to_sections[i] #sections associated with this element
    local ns = length(sections)     # number of sections in this element 

    #start from the middle-ish section
    mid = div(ns,2) #note that section is in the form of 1 to ns, consecutively.

    #get the feasible designs for the middle section
    feasible_idx = output_results[sections[mid]]
    
    #catalog was already sorted, so I think we can leave this part, just filter, to save time.
    sub_catalog = sort(catalog[feasible_idx, :], [:carbon,:fc′, :ec])
    # sub_catalog = catalog[feasible_idx, :]

    #now, loop each design in the sub catalog, see if as and fpe are available in all sections.
    #if not, remove that design from the sub catalog.
    #if yes, keep it.

    for d in eachrow(sub_catalog)
        # all_fc′ = true
        all_as = true
        all_fpe = true
        for s in sections
            #check if the design is available in that section.
            #if not, remove it from the sub catalog.
            #if yes, keep it.
            # if !(d[:fc′] ∈ catalog[output_results[s], :fc′])
            #     all_fc′ = false
            # end
            if !(d[:as] ∈ catalog[output_results[s], :as])
                all_as = false
            end
            if !(d[:fpe] ∈ catalog[output_results[s], :fpe])
                all_fpe = false
            end
        end
        if all_as && all_fpe
            #first run, found, move on.
            break
        end
    end

    # sort!(sub_catalog, [:carbon, :fpe, :as, :ec])

    #get the first one, they will appear in the entire thing anyway.
    this_fpe = sub_catalog[1,:fpe]
    this_as = sub_catalog[1, :as]

    sections_designs = Vector{Vector}(undef, ns)
    for is in eachindex(elements_to_sections[i])
        #current section index
        s = elements_to_sections[i][is]

        feasible_idx = output_results[s]
        sub_catalog = catalog[feasible_idx, :]

        fpe_as(fpe::Float64, as::Float64) = fpe == this_fpe && as == this_as

        this_catalog = filter([:fpe, :as] => fpe_as , catalog[output_results[s],:])

        sort!(this_catalog, [:carbon,:ec])


        
        #get the first one, it's the best.
        select_ID= this_catalog[1,:ID]
        #find lowest e for this one.
        sections_designs[is] = collect(catalog[select_ID,:])
        println(sections_designs[is])
        end

    elements_designs[i] = sections_designs

end

return elements_designs
end
        
elements_designs = find_optimum(output_results, demands)

list_fc′ = Vector{Float64}(undef, ns)
for i in keys(elements_designs)
    println(" #### ",i)
    for j in eachindex(elements_designs[i])
        println(elements_designs[i][j][9])
        println(elements_designs[i][j][1])
        # list_fc′[elements_designs[i][j][9]] = elements_designs[i][j][1]
    end
end


open("designs.json","w") do f
    JSON.print(f, elements_designs)
end