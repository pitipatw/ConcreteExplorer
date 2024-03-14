μs = 0.3
ne = length(elements_designs)
# for e in 1:ne
e = 1
    element_designs = elements_designs[e]
    ns = length(element_designs)
    support_dps = element_designs[1][:dps]
    next_dps = element_designs[2][:dps]

    another_support_dps = element_designs[ns][:dps]
    another_next_dps = element_designs[ns-1][:dps]

    as = element_designs[1][:as]
    fpe = element_designs[1][:fpe]
    

    
    #There is a symmetry (there must be a symmetry)
    @assert support_dps == another_support_dps "Symmetry Error at supports: $support_dps ≠ $another_support_dps."
    @assert next_dps == another_next_dps "Symmetry Error next to supports: $next_dps ≠ $another_next_dps."
    
    L = 500.0 # That's the distance between the 2 deviators.
    θ = atan((next_dps-support_dps)/L)
    rad2deg(θ)
    axial_component = as*fpe*cos(θ)
    maxV = 0
    for i in eachindex(element_designs)
        println(element_designs[i][:Vu])
        if element_designs[i][:Vu] > maxV 
            maxV = element_designs[i][:Vu]
        end
    end
    
    required_normal_force = maxV/μs
    axial_component = as*fpe*cos(θ)
    additional_force = required_normal_force - axial_component

    for s in eachindex(element_designs)
        element_designs[s][:axial_force] = additional_force
    end


