abstract type AbstractStructuralSection end 
abstract type AbstractPixelFrameSection <: AbstractStructuralSection end
abstract type AbstractReinforcedConcreteSection <: AbstractStructuralSection end 

mutable struct PixelFrameSection <: AbstractPixelFrameSection
    
    material::Material
    properties::PixelFrameProperties
    
    section::CompoundSection
    fc′::Float64
    fR1::Float64
    fR3::Float64

    pt_area::Vector{Float64}
    pt_force::Vector{Float64}
    pt_pos::Matrix{Float64}

    #Material Properties
    Ec::Float64  
    Ept::Float64 #Modulus of the Post Tension steel

    function PixelFrameSection(L::Real, t::Real, Lc::Real, fc′::Real,fR1::Float64, fR3::Float64, pt_area::Vector{Float64}, pt_force::Vector{Float64}, pt_pos::Matrix{Float64})
        @assert size(pt_pos)[1] == 2 "Post tension steel position must be a 2xn Matrix{Float64}"
        @assert length(pt_area) === length(pt_force) "Mis-match numbers of post tensioned areas and forces"
        @assert length(pt_area) === size(pt_pos)[2] "Mis-match numbers of post tensioned areas and positions"
        @assert length(pt_force) === size(pt_pos)[2] "Mis-match numbers of post tensioned forces and positions"
        
        section = make_Y_layup_section(L,t,Lc)
        Ec = 4700*sqrt(fc′) #MPa
        Ept = 200_000.0 #MPa
        pixelframesection = new(section,fc′, fR1, fR3, pt_area, pt_force, pt_pos,Ec, Ept)
        return pixelframesection
    end
end 

function Base.show(io::IO, p::PixelFrameSection)
    text = "PixelFrameSection\n"
    # fields = fieldnames(::PixelFrameSection)
    # Have to find a way to correctly/automatically match units with fields

    fields = (:section, :fc′, :fR1, :fR3, :pt_area, :pt_force, :pt_pos, :Ec, :Ept)
    units  = ("-", "MPa", "MPa", "MPa", "mm²", "N", "mm", "MPa", "MPa")
    for i  in eachindex(fields) 
        sub_text = string("    ", fields[i], ": ", getfield(p, fields[i])," ",units[i], "\n")
        text = text*sub_text 
    end 

    text = replace(text, ";" => "\n"*repeat(" ", 12))
    print(io,text)
end
