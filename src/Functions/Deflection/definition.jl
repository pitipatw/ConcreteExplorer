#struct definition
using AsapToolkit
# Constructing types
# Material Properties
mutable struct Material
    fc′::Float64 # Concrete strength [MPa] ****Should update on the test day using cylinder test***
    Ec::Float64 # MPa  ACI fc-> Concrete modulus relationship [MPa]
    Eps::Float64 #Post tensioning steel modulus [MPa]
    fpy::Float64 #MPa  
    #Safe load on the website https://www.engineeringtoolbox.com/wire-rope-strength-d_1518.html 
    # is ~ 150 MPa. Currently 140 MPa :)
end

mutable struct Section
    em::Float64 # Eccentricity at the middle of the member [mm]
    es::Float64 # Eccentricity at the support of the member[mm]
    em0::Float64 # Initial eccentricity at the midspan [mm]
    dps0::Float64 # Initial distance from the top to the point of application of the load [mm]
    Ls::Float64 # Distance from support to the first load point [mm]
    Ld::Float64 # Distance from support to the first deviator [mm]
    L::Float64  # Total length of the member [mm]
    # two 1/4" bars with 1200 lb capacity
    Aps::Float64 # Total area of the steel in the section [mm^2]
    Atr::Float64 # Transformed area of the cross section [mm^2]
    Itr::Float64 # Moment of inertia of the transformed cross section [mm^4]
    Zb::Float64  # Section modulus of the concrete section from the centroid to extreme tension fiber [mm^3]
    K1::Float64
    K2::Float64

    function Section(em::F)

        if Ld < Ls
            K1 = Ls / L - 1
            K2 = Ls / Ls(Ld / L)^2 - (Ls / L)^2
        elseif Ld ≥ Ls
            K1 = Ld / L - 1
            K2 = 0
        end
    end

end

mutable struct Loads
    w::Float64 # Selfweight [N/mm]
    Mg::Float64 # Moment due to selfweight [Nmm]
    fr::Float64 # Concrete cracking strenght [MPa]
    r::Float64 # Radius of gyration [mm]
    #ps_force::Float64 # Post tensioning force [N]
    fpe::Float64 # Effective post tensioning stress [MPa]
    ϵpe::Float64 # Effective post tensioning strain [-]
    ϵce::Float64 # Effective concrete strain [-]
    Mdec::Float64 # decompression moment [Nmm]
end

mutable struct Element
    Mat::Material
    Sec::Section
    F::Loads
    geometry::CompoundSection
    type::String
end
