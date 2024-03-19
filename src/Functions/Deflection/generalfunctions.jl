include("../Geometry/pixelgeo.jl")

abstract type AbstractStructuralElement end 
abstract type AbstractPixelFrameElement <: AbstractStructuralElement end 

abstract type AbstractStructuralSection end 
abstract type AbstractPixelFrameSection <: AbstractStructuralSection end

mutable struct PixelFrameSection <: AbstractPixelFrameSection
    section::CompoundSection
    fc′::Float64
    fR1::Float64
    fR3::Float64

    pt_area::Vector{Float64}
    pt_force::Vector{Float64}
    pt_pos::Matrix{  }

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

"""
    PixelFrameElement
Create an element consists of PixelFrame sections
"""
mutable struct PixelFrameElement <: AbstractPixelFrameElement
    # sections::Vector{AbstractPixelFrameSection}
    sections::Vector{CompoundSection}
    fc′::Float64 # Concrete strength [MPa] ****Should update on the test day using cylinder test***
    Ec::Float64 # MPa  ACI fc-> Concrete modulus relationship [MPa]
    Eps::Float64 #Post tensioning steel modulus [MPa]
    fpy::Float64 #MPa  

    em::Float64 # Eccentricity at the middle of the member [mm]
    es::Float64 # Eccentricity at the support of the member   [mm]
    em0::Float64 # Initial eccentricity at the midspan        [mm]
    dps0::Float64 # Initial distance from the top to the point of application of the load [mm]
    Ls::Float64 # Distance from support to the first load point [mm]
    Ld::Float64 # Distance from support to the first deviator [mm]
    L::Float64 # Total length of the member [mm]
    # two 1/4" bars with 1200 lb capacity
    Aps::Float64 # Total area of the steel in the section [mm^2]
    Atr::Float64 # Transformed area of the cross section <at the mid span?> [mm^2]
    Itr::Float64 # Moment of inertia of the transformed cross section [mm^4]
    Zb::Float64 # Section modulus of the concrete section from the centroid to extreme tension fiber [mm^3]


    w::Float64 # Selfweight [N/mm]
    mg::Float64 # Moment due to selfweight [Nmm]
    fr::Float64 # Concrete cracking strenght [MPa]
    r::Float64 # Radius of gyration [mm]
    ps_force::Float64 # Post tensioning force [N]
    concrete_force::Float64
    fpe::Float64 # Effective post tensioning stress [MPa]
    ϵpe::Float64 # Effective post tensioning strain [-]
    ϵce::Float64 # Effective concrete strain [-]
    Mdec::Float64 # decompression moment [Nmm]

    test::Bool#True if this is from the half-scale test.

        """PixelFrameElement(L::Float64, t::Float64, Lc::Float64, 
        fc′::Float64, as::Float64, fpe::Float64, Le::Float64, dps0::Float64,) 
        """


        """ PixelFrameElement()
        Input cases for the half-scale test
        """
    function PixelFrameElement()
        println("Creating a PixelFrame element from the half-scale test")
        println("Notes that numbers from the field 'compoundsection' might be different from its actual properties")
        println("This is due to the holes and cuts for the actual test")

        compoundsection = make_Y_layup_section(141.7,19.1,25.2)
        pixelframeelement = new([compoundsection]) #it takes a vector of compound section

        #..........Notes..........
        # Use Ld = Ls (this test only) 
        # Eccentricities measured from the neutral axis
        # M is the moment in the constant region
        # Mg = moment due to the selfweight
        # M(x) is the moment equation due to the load
        #Units N, mm, MPa

        #   Material Properties
        pixelframeelement.fc′ = 36.0 # Concrete strength [MPa] ****Should update on the test day using cylinder test***
        # Ec = 4700.0*sqrt(fc′) # MPa  ACI fc-> Concrete modulus relationship [MPa]
        pixelframeelement.Ec = 58000.0 # MPa  from the cylinder test
        pixelframeelement.Eps = 70000.0 #Post tensioning steel modulus [MPa]
        pixelframeelement.fpy = 0.002 * pixelframeelement.Eps #MPa  
        #Safe load on the website https://www.engineeringtoolbox.com/wire-rope-strength-d_1518.html 
        # is ~ 150 MPa. Currently 140 MPa :)

        # PixelFrame section/element properties
        centroid_to_top = 91.5 #[mm]
        pixelframeelement.em = 230.0 # Eccentricity at the middle of the member [mm]
        pixelframeelement.es = 0.0 # Eccentricity at the support of the member   [mm]
        pixelframeelement.em0 = 230.0 # Initial eccentricity at the midspan        [mm]

        pixelframeelement.dps0 = centroid_to_top + pixelframeelement.em0 # Initial distance from the top to the point of application of the load [mm]
        pixelframeelement.Ls = 502.7 # Distance from support to the first load point [mm]
        pixelframeelement.Ld = 502.7 # Distance from support to the first deviator [mm]
        pixelframeelement.L = 2000.0 # Total length of the member [mm]
        # two 1/4" bars with 1200 lb capacity
        pixelframeelement.Aps = 2.0 * (0.25 * 25.4)^2 * pi / 4.0 # Total area of the post tensioned steel [mm2]
        #Pure concrete area = 18537.69 mm2
        #Transformed steel area = 347.96 mm2 
        pixelframeelement.Atr = 18537.69 # Transformed area of the cross section [mm2] (= Concrete area if there is no embedded rebars)
        pixelframeelement.Itr = 6.4198e+07 #moment of inertia [mm4], no embedded steel, therefore, only from concrete.
        # pixelframeelement.Itr = 1.082e+8 #this number includes deviated steels.

        pixelframeelement.Zb = pixelframeelement.Itr/centroid_to_top # Elastic modulus of the concrete section from the centroid to extreme tension fiber [mm3]
        # If there are multiple materials, transformed section geometry is needed for Zb (and everything related to section area)

        #forces
        pixelframeelement.w = pixelframeelement.Atr / 10^9 * 2400.0 * 9.81 # Selfweight [N/mm]
        pixelframeelement.mg = pixelframeelement.w * pixelframeelement.L^2 / 8.0 # Moment due to selfweight [Nmm]
        pixelframeelement.fr = 0.7 * sqrt(pixelframeelement.fc′) # Concrete cracking strenght [MPa]
        pixelframeelement.r = sqrt(pixelframeelement.Itr / pixelframeelement.Atr) # Radius of gyration [mm]
        pixelframeelement.ps_force = 890.0/sind(24.0) # Post tensioning force [N]
        pixelframeelement.Mdec = pixelframeelement.ps_force*pixelframeelement.em
        pixelframeelement.concrete_force = pixelframeelement.ps_force*cos(24.0*pi/180.0) # 
        pixelframeelement.fpe = pixelframeelement.ps_force/pixelframeelement.Aps # Effective post tensioning stress [MPa] ***will input the one on the test day***
        pixelframeelement.ϵpe = pixelframeelement.fpe / pixelframeelement.Eps # Effective post tensioning strain [mm/mm]
        #find moment due to the applied force.
        pixelframeelement.ϵce = pixelframeelement.ps_force*pixelframeelement.em/pixelframeelement.Zb/pixelframeelement.Ec - pixelframeelement.concrete_force/pixelframeelement.Atr/pixelframeelement.Ec # effetive strain in the concrete [mm/mm]
        #for using test setup
        pixelframeelement.test = true
        return pixelframeelement
    end

    """
        (To do)Create PixelFrame element from a vector of PixelFrameSection
    """
    function pixelframeelement(pixelframesections::Vector{PixelFrameSection})
        sections = getfield.(pixelframesections,section)

        compoundsection = sections[1]
        pixelframeelement = new([compoundsection]) #it takes a vector of compound section

        #..........Notes..........
        # Use Ld = Ls (this test only) 
        # Eccentricities measured from the neutral axis
        # M is the moment in the constant region
        # Mg = moment due to the selfweight
        # M(x) is the moment equation due to the load
        #Units N, mm, MPa

        #   Material Properties
        pixelframeelement.fc′ = 36.0 # Concrete strength [MPa] ****Should update on the test day using cylinder test***
        # Ec = 4700.0*sqrt(fc′) # MPa  ACI fc-> Concrete modulus relationship [MPa]
        pixelframeelement.Ec = 58000.0 # MPa  from the cylinder test
        pixelframeelement.Eps = 70000.0 #Post tensioning steel modulus [MPa]
        pixelframeelement.fpy = 0.002 * pixelframeelement.Eps #MPa  
        #Safe load on the website https://www.engineeringtoolbox.com/wire-rope-strength-d_1518.html 
        # is ~ 150 MPa. Currently 140 MPa :)

        # PixelFrame section/element properties
        centroid_to_top = 91.5 #[mm]
        pixelframeelement.em = 230.0 # Eccentricity at the middle of the member [mm]
        pixelframeelement.es = 0.0 # Eccentricity at the support of the member   [mm]
        pixelframeelement.em0 = 230.0 # Initial eccentricity at the midspan        [mm]

        pixelframeelement.dps0 = centroid_to_top + pixelframeelement.em0 # Initial distance from the top to the point of application of the load [mm]
        pixelframeelement.Ls = 502.7 # Distance from support to the first load point [mm]
        pixelframeelement.Ld = 502.7 # Distance from support to the first deviator [mm]
        pixelframeelement.L = 2000.0 # Total length of the member [mm]
        # two 1/4" bars with 1200 lb capacity
        pixelframeelement.Aps = 2.0 * (0.25 * 25.4)^2 * pi / 4.0 # Total area of the post tensioned steel [mm2]
        #Pure concrete area = 18537.69 mm2
        #Transformed steel area = 347.96 mm2 
        pixelframeelement.Atr = 18537.69 # Transformed area of the cross section [mm2] (= Concrete area if there is no embedded rebars)
        pixelframeelement.Itr = 6.4198e+07 #moment of inertia [mm4], no embedded steel, therefore, only from concrete.
        # pixelframeelement.Itr = 1.082e+8 #this number includes deviated steels.

        pixelframeelement.Zb = pixelframeelement.Itr/centroid_to_top # Elastic modulus of the concrete section from the centroid to extreme tension fiber [mm3]
        # If there are multiple materials, transformed section geometry is needed for Zb (and everything related to section area)

        #forces
        pixelframeelement.w = pixelframeelement.Atr / 10^9 * 2400.0 * 9.81 # Selfweight [N/mm]
        pixelframeelement.mg = pixelframeelement.w * pixelframeelement.L^2 / 8.0 # Moment due to selfweight [Nmm]
        pixelframeelement.fr = 0.7 * sqrt(pixelframeelement.fc′) # Concrete cracking strenght [MPa]
        pixelframeelement.r = sqrt(pixelframeelement.Itr / pixelframeelement.Atr) # Radius of gyration [mm]
        pixelframeelement.ps_force = 890.0/sind(24.0) # Post tensioning force [N]
        pixelframeelement.Mdec = pixelframeelement.ps_force*pixelframeelement.em
        pixelframeelement.concrete_force = pixelframeelement.ps_force*cos(24.0*pi/180.0) # 
        pixelframeelement.fpe = pixelframeelement.ps_force/pixelframeelement.Aps # Effective post tensioning stress [MPa] ***will input the one on the test day***
        pixelframeelement.ϵpe = pixelframeelement.fpe / pixelframeelement.Eps # Effective post tensioning strain [mm/mm]
        #find moment due to the applied force.
        pixelframeelement.ϵce = pixelframeelement.ps_force*pixelframeelement.em/pixelframeelement.Zb/pixelframeelement.Ec - pixelframeelement.concrete_force/pixelframeelement.Atr/pixelframeelement.Ec # effetive strain in the concrete [mm/mm]
        #for using test setup
        pixelframeelement.test = true
        return pixelframeelement
    end
end





"""
Currently unused.
"""
function pixelframe_properties!(pixelframesection::AbstractPixelFrameSection)::AbstractPixelFrameSection

    pixelframesection.fc′::Float64 = 0.0
    pixelframesection.fR1::Float64 = 0.0
    pixelframesection.fR3::Float64 = 0.0

    pixelframesection.pt_area::Vector{Float64} = [0.0]
    pixelframesection.pt_force::Vector{Float64} = [0.0]
    pixelframesection.pt_pos::Matrix{Float64} = [0.0 0.0]

    #Material Properties
    pixelframesection.Ec::Float64  = 0.0
    pixelframesection.Ept::Float64 = 0.0


    return pixelframesection
end


