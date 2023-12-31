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

"""

"""
function beam_design(fc′::Float64, ec_concrete::Float64;
    span::Float64 = 10_000.0,
    bay::Float64  = 10_000.0,
    # wLL::Float64  = 1.92e-3,
    wLL::Float64  = 2.40e-3,
    wSDL::Float64 = (24e-6*300), #2400kg/m3 * 300cm -> N/mm2
    echo = false
    )

    ρ_concrete = 24000.0e-9 #N/mm3
    Ec = 4700*sqrt(fc′) #MPa
    fy = 420.0 #MPa
    Es = 210_000 #MPa
    """
    Define a group of rebars for a RC section
    ecc of rebar in from CLF database -> fabricated rebars
    (unfabricated)
    753 kgCO2e/metricTon -> 753*7.85 metricTon/m3 = 5911.05 kgCO2e/m3
    or
    (fabricated)
    854 kgCO2e/metricTon -> 854*7.85 metricTon/m3 = 6703.90 kgCO2e/m3
    """
    ec_rebar = 5911.05e-9 #kgCO2e/mm3
    #section definition
    d = 600.0 #mm
    bd_ratio = 0.5
    b = bd_ratio*d #mm
    covering = 50.0 #mm 
    p1 = [0.0 , 0.0]
    p2 = [b, 0.0]
    p3 = [b, -d-covering]
    p4 = [0.0, -d-covering]

    points = [p1, p2, p3, p4]
    section = SolidSection(points)

    # d_cal = section.ymax -section.ymin
    # @assert d_cal == d
    #have to find the least amount of embodied carbon for the

    w_load = (1.2*wSDL + 1.6*wLL)*bay
    w_beamDL = 1.2*ρ_concrete*section.area
    w = w_load + w_beamDL
    Mdemand = w * span^2/8
    Vdemand = w * span/2
    # println("Moment Demand: ", Mdemand)
    # println("Shear Demand: ", Vdemand)
    try
        C = 2*Mdemand/(0.85*fc′*b)
        B = -2*d
        a = (-B - sqrt(B^2-4*C))/2
    catch 
        return nothing, false
    end

    if a < 0 
        # println("a is less than 0")
        a = (-B + sqrt(B^2-4*C))/2
    end

    @assert 0 < a < d
    @assert  Mdemand - 0.85*fc′*a*b*(d-a/2) <= 0.001


    as = Mdemand/(fy*(d-a/2))
    c = a/0.85
    ϵs = 0.003*(d-c)/c 
    if ϵs < 0.005
        # println("NON ductile behavior")
    end
    #0.0018 is from ACI (see the jupyter notebook)
    ρ_min = clamp(as/(section.area-as), 0.0018, Inf)

    #find ρ_max 
    c_max = 3/8*d
    a_max = 0.85*c_max
    as = Mdemand/(fy*(d-a_max/2))
    ρ_max = as/(section.area-as)

    #here, should we select the ρ_min always? 
    ρ_selected = (ρ_min+ρ_max)/2
    area_max = section.area*ρ_max
    area_min = section.area*ρ_min
    
    a_maximum = area_max*fy/(0.85*fc′*b)
    c_maximum = a_maximum/0.85
    ϵs_maximum = 0.003*(d-c_maximum)/c_maximum

    a_minimum = area_min*fy/(0.85*fc′*b)
    c_minimum = a_minimum/0.85
    ϵs_minimum = 0.003*(d-c_minimum)/c_minimum


    reinforcement_area = section.area*ρ_selected
    a = reinforcement_area*fy/(0.85*fc′*b)
    c = a/0.85
    ϵs = 0.003*(d-c)/c
    if ϵs < 0.005
        # println("Mean reinforcement is not ductile enough, using ρ_minimum instead")
        ρ_selected = ρ_min
        reinforcement_area = section.area*ρ_selected
        a = reinforcement_area*fy/(0.85*fc′*b)
        c = a/0.85
        ϵs = 0.003*(d-c)/c
    end

    # δ_max = 5*w*L^4/(384*E*I)

    #find E and I of the beam. 
    #we turn everything into a pure concrete with E = 4700sqrt(fc′) 
    as_transformed = Es/Ec*as 
    total_area = section.area-as+as_transformed
    #find a new neutral axis of the section 
    NA_concrete = section.centroid[2]
    NA_rebar = -d
    NA_section = (section.area*NA_concrete+as_transformed*NA_rebar)/total_area
    @assert -d < NA_section < NA_concrete

    I_RCsection = section.Ix + section.area*(NA_section-NA_concrete)^2 + as_transformed*(NA_rebar-NA_section)^2

    @assert section.Ix < I_RCsection

    δ_max = 5*w*span^4/(384*Ec*I_RCsection)
    #Deflection limit L/240 
    δ_limit =  span/240

    # @assert δ_max <= δ_limit
    serviceability = false
    if δ_max <= δ_limit
        serviceability = true
    end

    #Flexural design is done. 
    if echo
        @printf "Concrete section %d x %d sqmm\n" section.xmax-section.xmin section.ymax-section.ymin
        @printf "Reinforcement area: %.2f sqmm\n" reinforcement_area
        println("Serviceability: ", serviceability ? "Pass!" : "Fail")
    end


    #create points on circle for 2 rebars
    offsety = -d 
    offsetx = b/4
    a = reinforcement_area/2 
    r = sqrt(a/pi)
    pts1 = [[r*cos(θ)+offsetx, r*sin(θ)+offsety] for θ in 0:2*pi/10:2pi]
    pts2 = [[r*cos(θ)+b-offsetx, r*sin(θ)+offsety] for θ in 0:2*pi/10:2pi]

    circle1 = SolidSection(pts1)
    circle2 = SolidSection(pts2)

    embodied_carbon = span*(ec_concrete*(section.area-reinforcement_area) + ec_rebar*reinforcement_area)
    rc_section = create_rc_section(section, [circle1, circle2], d, embodied_carbon)
    return rc_section, serviceability
end

