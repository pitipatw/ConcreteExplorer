using AsapSections

mutable struct PixelFrameSection <: StructuralSection
    section::Section
    fc′::Float64
    fR1::Float64
    fR3::Float64
    as::Vector{Float64}
    fps::Vector{Float64}

    
end 

mutable struct ReinforcedConcreteSection <: StructuralSection
end 

