using Symbolics
using Roots

"""
Tester
find_fc′ should ≈ fc′00
"""
# b1 = 200.0
# d1 = 200.0
# Mu1 = 2.0e6
# As1 = 800.0
# fy1 = 420.0
# function getmu(fc′1)
#     a = As1 * fy1 / (0.85 * fc′1 * b1)
#     β1 = clamp(0.85 - 0.05 * (fc′1 - 28) / 7, 0.65, 0.85)
#     c = a / β1
#     ϵs = (d1 - c) / c * 0.003
#     ϕ = clamp(0.65 + 0.25 * (ϵs - 0.002) / 0.003, 0.65, 0.90)
#     return Mu1 = ϕ * As1 * fy1 * (d1 - a / 2)
# end

# fc′00 = 35
# Mu0 = getmu(fc′00)

function define_function(As0::Float64, fy0::Float64, b0::Float64, d0::Float64, Mu0::Float64)
    @variables Mu, As, fy, b, d, fc′
    a = As * fy / (0.85 * fc′ * b)
    β1 = clamp(0.85 - 0.05 * (fc′ - 28) / 7, 0.65, 0.85)
    c = a / β1
    ϵs = 0.003 * (d - c) / c
    ϕ = clamp(0.65 + 0.25 * (ϵs - 0.002) / 0.003, 0.65, 0.90)
    eq = Mu - ϕ * As * fy * (d - a / 2)
    simp_eq = simplify(eq)
    # @show simp_eq
    subs = Dict(Mu => Mu0, As => As0, fy => fy0, b => b0, d => d0)
    f22 = substitute(simp_eq, subs)
    ff = eval(build_function(f22, fc′))
    return ff
end

# ff = define_function(As1, fy1, b1, d1, Mu0)

function find_fc′(ff; fc′0 = 28.0)
    return find_zero(ff, fc′0)
end

# find_fc′(ff) 

# function find_fc′(As0::Float64, fy0::Float64, b0::Float64, d0::Float64, Mu0::Float64; fc′0 = 28.0)
#     ff = define_function(As1, fy1, b1, d1, Mu0)
#     return find_zero(invokelatest(ff,0), fc′0)
# end

# find_fc′(As1, fy1, b1, d1, Mu0)