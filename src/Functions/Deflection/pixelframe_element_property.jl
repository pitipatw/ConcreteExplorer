abstract type PixelFrameElementProperty end
abstract type AbstractPixelFrameElementProperty <: PixelFrameElementProperty end 


"""
function Material(fc′::Real, Ec::Real, Eps::Real, fpy::Real, fr::Real)::Material <: AbstractPixelFrameElementProperty\n
function Material(fc′::Real, Eps::Real, fpy::Real)::Material <: AbstractPixelFrameElementProperty\n
    Ec::Float64 = 4700*sqrt(fc′) #MPa\n
    fr::Float64 = 0.62*sqrt(fc′) #MPa\n
\n
function Material(fc′::Real)::Material <: AbstractPixelFrameElementProperty\n
    Ec::Float64  = 4700*sqrt(fc′) #MPa\n
    Eps::Float64 = 200_000        #MPa\n
    fpy::Float64 = 420            #MPa\n
    fr::Float64  = 0.62*sqrt(fc′) #MPa\n
"""
mutable struct Material <: AbstractPixelFrameElementProperty
    fc′::Float64    # Concrete strength [MPa]
    Ec::Float64     # MPa  ACI fc-> Concrete modulus relationship [MPa]
    Eps::Float64    # Post tensioning steel modulus [MPa]
    fpy::Float64    # yield strength of rebar [MPa]
    fr::Float64     # Concrete cracking strenght [MPa]

    """
    function Material(fc′::Real, Ec::Real, Eps::Real, fpy::Real, fr::Real)::Material <: AbstractPixelFrameElementProperty
    """
    function Material(fc′::Real, Ec::Real, Eps::Real, fpy::Real, fr::Real)
        return new(fc′, Ec, Eps, fpy,fr)
    end

    """
    function Material(fc′::Float64, Eps::Float64, fpy::Float64)::Material <: AbstractPixelFrameElementProperty
        Ec::Float64 = 4700*sqrt(fc′) #MPa
        fr::Float64 = 0.62*sqrt(fc′) #MPa
    """
    function Material(fc′::Real, Eps::Real, fpy::Real)
        Ec::Float64 = 4700*sqrt(fc′) #MPa
        fr::Float64 = 0.62*sqrt(fc′) #MPa

        return new(fc′, Ec, Eps, fpy, fr)
    end

    """
    function Material(fc′:Real)::Material <: AbstractPixelFrameElementProperty
        Ec::Float64  = 4700*sqrt(fc′) #MPa
        Eps::Float64 = 200_000        #MPa
        fpy::Float64 = 420            #MPa
        fr::Float64  = 0.62*sqrt(fc′) #MPa
    """
    function Material(fc′::Real)
        Ec::Float64  = 4700*sqrt(fc′) #MPa
        Eps::Float64 = 200_000        #MPa
        fpy::Float64 = 420            #MPa
        fr::Float64  = 0.62*sqrt(fc′) #MPa

        return new(fc′, Ec, Eps, fpy, fr)
    end
end

mutable struct Properties <: AbstractPixelFrameElementProperty #PixelFrameProperties
    em::Float64     # Eccentricity at the middle of the member [mm]
    es::Float64     # Eccentricity at the support of the member[mm]
    θ::Float64     # Angle of the tendon at the support, measured from the element's horizontal line
    em0::Float64    # Initial eccentricity at the midspan [mm]
    dps0::Float64   # Initial distance from the top to the point of application of the load [mm]
    L::Float64      # Total length of the member [mm]
    Ls::Float64     # Distance from support to the first load point [mm]
    Ld::Float64     # Distance from support to the first deviator [mm]
    Aps::Float64    # Total area of the steel in the section [mm^2]
    Atr::Float64    # Transformed area of the cross section [mm^2]
    Itr::Float64    # Moment of inertia of the transformed cross section [mm^4]
    r::Float64      # Radius of gyration [mm]
    Zb::Float64     # Section modulus of the concrete section from the centroid to extreme tension fiber [mm^3]
    
    w::Float64                    # Selfweight [N/mm]
    moment_selfweight::Float64    # Moment due to selfweight [Nmm]
    fpe::Float64                  # Effective post tensioning stress [MPa]
    moment_decompression::Float64 # decompression moment [Nmm]
    ϵpe::Float64                  # Effective post tensioning strain [-]
    ϵce::Float64                  # Effective concrete strain due to the post tension force[-]

    # K1::Float64     # A parameter from Ng and Tan paper [-]
    # K2::Float64     # A parameter from Ng and Tan paper[-]

    function Properties(compoundsection::CompoundSection, material::Material, Aps::Real, θ::Real, L::Real, Ls::Real, Ld::Real,em::Float64, es::Float64, fpe::Real)
        @unpack Ec, Eps = material
        em0::Float64  = em 
        dps0::Float64 = em
        Atr::Float64  = compoundsection.area #Transformed section between concrete and "nonprestressed steel".
        Itr::Float64  = compoundsection.Ix
        r::Float64    = sqrt(Itr/Atr)
        Zb::Float64   = Itr/(compoundsection.centroid[2] - compoundsection.ymin)
        
        w::Float64    = 2400*9.81*compoundsection.area/1e9 #N/mm
        moment_selfweight::Float64 = w*L^2/8
        ps_force = Aps*fpe
        moment_decompression::Float64 = ps_force*em0
        ϵpe::Float64 = ps_force*em0/Eps
        ϵce::Float64 = (ps_force*em/Zb - ps_force*cos(θ)/Atr)/Ec

        # if Ld < Ls
        #     K1 = Ls / L - 1
        #     K2 = Ls / Ls(Ld / L)^2 - (Ls / L)^2
        # elseif Ld ≥ Ls
        #     K1 = Ld / L - 1
        #     K2 = 0
        # end

        # @assert L == Ls + Ld + Ls "Error, total length (L) ≠ 2Ls + Ld"
        # return new(em,es,em0,dps0,Ls,Ld,L, Aps, Atr, Itr, Zb, K1, K2) 

        return new(em, es, θ, em0, dps0, 
            L, Ls, Ld, Aps, Atr, Itr, r, Zb,
            w,moment_selfweight, fpe,moment_decompression, 
            ϵpe, ϵce)
    end
end


