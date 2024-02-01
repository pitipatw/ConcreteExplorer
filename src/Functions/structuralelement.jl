using AsapSections

include("Geometry/pixelgeo.jl")
include("GeneralFunctions.jl")

#Structural Element
abstract type StructuralElement end 
abstract type ReinforcedConcreteElement <: StructuralElement  end 
abstract type PixelFrameElement <: StructuralElement end 



#Structural Section
abstract type StructuralSection end 
abstract type ReinforcedConcreteSection <: StructuralSection end 
abstract type PixelFrameSection <: StructuralSection end 

mutable struct PixelFrameSection
    section::CompoundSection
    config::String

    fc′::Float64
    fR1::Float64
    fR3::Float64

    pt_area::Vector{Float64}
    pt_force::Vector{Float64}
    pt_pos::Matrix{Float64}

    #Material Properties
    Ec::Float64 
    Eps::Float64

    """
        PixelFrameSection(L::Float64, t::Float64, Lc::Float64; config::String = "Y")
    Create a PixelFrame section from L, t, and Lc parameters.
    """
    function PixelFrameSection(L::Float64, t::Float64, Lc::Float64; config::String = "Y")
        compoundsection = make_Y_layup_section(L,t,Lc)
        pixelframesection = new(compoundsection, config)
        
        #populate PixelFrameSection properties
        pixelframe_properties!(pixelframesection)

        return pixelframesection
    end
end 


"""
    PixelFrameElement
Create an element consists of PixelFrame sections
"""
mutable struct PixelFrameElement
    sections::Vector{PixelFrameSection}

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
    #ps_force::Float64 # Post tensioning force [N]
    fpe::Float64 # Effective post tensioning stress [MPa]
    ϵpe::Float64 # Effective post tensioning strain [-]
    ϵce::Float64 # Effective concrete strain [-]
    Mdec::Float64 # decompression moment [Nmm]

        """ PixelFrameElement(test::Bool)
        Input cases for the half-scale test
        """
    function PixelFrameElement(test::Bool)
        compoundsection = 0.0
        pixelframeelement = new([compundsection]) #it takes a vector of compound section

        #..........Notes..........
        # Use Ld = Ls (this test only) 
        # Eccentricities measured from the neutral axis
        # M is the moment in the constant region
        # Mg = moment due to the selfweight
        # M(x) is the moment equation due to the load
        #Units N, mm, MPa

        #   Material Properties
        pixelframelement.fc′ = 36.0 # Concrete strength [MPa] ****Should update on the test day using cylinder test***
        # Ec = 4700.0*sqrt(fc′) # MPa  ACI fc-> Concrete modulus relationship [MPa]
        pixelframelement.Ec = 58000.0 # MPa  from the cylinder test
        pixelframelement.Eps = 70000.0 #Post tensioning steel modulus [MPa]
        pixelframelement.fpy = 0.002 * Eps #MPa  
        #Safe load on the website https://www.engineeringtoolbox.com/wire-rope-strength-d_1518.html 
        # is ~ 150 MPa. Currently 140 MPa :)

        # PixelFrame section/element properties
        centroid_to_top = 91.5 #[mm]
        pixelframelement.em = 230.0 # Eccentricity at the middle of the member [mm]
        pixelframelement.es = 0.0 # Eccentricity at the support of the member   [mm]
        pixelframelement.em0 = em # Initial eccentricity at the midspan        [mm]
        pixelframelement.dps0 = centroid_to_top + em0 # Initial distance from the top to the point of application of the load [mm]
        pixelframelement.Ls = 502.7 # Distance from support to the first load point [mm]
        pixelframelement.Ld = Ls    # Distance from support to the first deviator [mm]
        pixelframelement.L = 2000.0 # Total length of the member [mm]
        # two 1/4" bars with 1200 lb capacity
        pixelframelement.Aps = 2.0 * (0.25 * 25.4)^2 * pi / 4.0 # Total area of the post tensioned steel [mm2]
    
        pixelframelement.Atr = 18537.69 # Transformed area of the cross section [mm2]
        pixelframelement.Itr = 6.4198e+07 #moment of inertia [mm4]
        # Itr = 1.082e+8 

        pixelframelement.Zb = Itr/centroid_to_top # Elastic modulus of the concrete section from the centroid to extreme tension fiber [mm3]
        # If there are multiple materials, transformed section geometry is needed for Zb (and everything related to section area)


        #forces
        pixelframelement.w = Atr / 10^9 * 2400.0 * 9.81 # Selfweight [N/mm]
        pixelframelement.mg = w * L^2 / 8.0 # Moment due to selfweight [Nmm]
        pixelframelement.fr = 0.7 * sqrt(fc′) # Concrete cracking strenght [MPa]
        pixelframelement.r = sqrt(Itr / Atr) # Radius of gyration [mm]
        @show pixelframelement.ps_force = 100 # 890.0/sind(24.0) # Post tensioning force [N]
        Mdec = ps_force*em
        @show concrete_force = ps_force*cos(24.0*pi/180.0) # 
        fpe = ps_force/Aps # Effective post tensioning stress [MPa] ***will input the one on the test day***
        ϵpe = fpe / Eps # Effective post tensioning strain [mm/mm]
        #find moment due to the applied force.
        ϵce = ps_force*em/Zb/Ec - concrete_force/Atr/Ec # effetive strain in the concrete [mm/mm]

        return
    end
end





mutable struct PixelFrameBeam 

end
