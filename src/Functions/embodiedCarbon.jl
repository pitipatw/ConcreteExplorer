#Functions associated with embodied carbon

"""
    fc2e(fc′::Real)::Float64 [kgCO2e/m3]
Return embodied carbon coefficient of concrete based on fc′ 


"""
function fc2e(fc′::Real; mode::String = "default")::Float64
    if mode =="default"
        #data on 25-80 MPa (Holcim) -> linearly interpolated to 100MPa.
        cfc = 1.6947*fc′ + 267.53
    elseif mode == "Holcim fiber"
        #only valid up to 60MPa.
        cfc = -0.0626944435544512 * fc′^2 + 10.0086510099949 * fc′ + 84.14807
    elseif mode == "myResults"
        #using data from J Broyles's paper
        cfc = 4.5726*fc′ + 217.29
    else 
        println("Invalid Mode")
        return nothing
    end

    return cfc
end


""" fc2e(fc′::Real, dosage::Real)::Float64 [kgCO2e/m3]
Get the embodied carbon coefficient of a concrete mix with steel fiber dosage
"""
function fc2e(fc′::Real, dosage::Real)::Float64
    cfc = fc2e(fc′, mode = "Holcim fiber") + 1.4*dosage #dosage: kg steel/m3 concret. 
    return cfc
end


"""
    rebar2e()::Float64
Get the embodied carbon coefficient for steel.
input : no input
output: ecc of rebar [kgCO2e/m3]
"""
function rebar2e()::Float64
    # 1.4 kgCO2e per kg steel (EC3)
    # 7850 kg/m3 ASTM A992 
    cst = 1.4*7850 #kgCO2e/m3
    return cst
end

using Makie

f_fc′_EC = Figure(size = (500,500)) 
ax1 = Axis(f_fc′_EC[1,1])

fc′ = 20:0.5:100
ec = fc2e.(fc′)

scatter!(ax1, fc′, ec)