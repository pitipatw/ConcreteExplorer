using AsapSections
using Makie, GLMakie
using Printf

mutable struct ReinforcedConcreteSection
    concrete_section::SolidSection
    rebar_section::Vector{SolidSection}
    depth::Float64
    embodied_carbon::Float64
    covering::Float64

    ec_concrete::Float64
    ec_rebar::Float64
end
function create_rc_section(
    concrete_section::SolidSection, 
    rebar_section::Vector{SolidSection},
    depth::Float64, 
    embodied_carbon::Float64)
    return ReinforcedConcreteSection(concrete_section, rebar_section,depth,embodied_carbon, 50.0,1.0,1.0)
end
function draw(rc_section::ReinforcedConcreteSection)
    f1 = Figure(size = (300,300)) 
    a1 = Axis(f1[1,1], aspect = DataAspect() )
    poly!(a1, rc_section.concrete_section.points, color = colorant"#B2B2B2")

    for i in rc_section.rebar_section
        poly!(a1, i.points, color = colorant"#3EA8DE")
    end

    return f1
end

function beam_design(Mdemand::Float64, fc′::Float64, ec::Float64)
    #have to find the least amount of embodied carbon for the



    return rc_section, serviceability
end

