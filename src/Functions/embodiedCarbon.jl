#Functions associated with embodied carbon

"""
Get embodied carbon coefficient of concrete based on fc′
input : fc′ [MPa]
output: ecc of fc′ [kgCO2e/m3]
"""
function fc2e(fc′::Real)
    out = -0.0626944435544512 * fc′^2 + 10.0086510099949 * fc′ + 84.14807
   return out
end

function rebar2e(area::Real)
    out = 120.0
    return out
end