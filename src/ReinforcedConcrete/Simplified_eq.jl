"""
Simplified equation (in psi and kips)
As = Mu/(4d)
"""
function simp_eq(Mu::Float64, d::Float64)
    As = Mu/4/d
    return As
end
