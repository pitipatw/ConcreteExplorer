using Nonconvex

function beam_obj(x)
    return x.^2 
end

function beam_cons(input_vector::Vector{Float64};
    fc′::Float64 = fc′, 
    ec_concrete::Float64 = ec,
    span::Float64 = 10_000.0,
    bay::Float64  = 10_000.0,
    # wLL::Float64  = 1.92e-3,
    wLL::Float64  = 2.40e-3,
    wSDL::Float64 = (300*24e-6), #2400kg/m3 * 300cm -> N/mm2
    ρ_concrete = 24000.0e-9, #N/mm3
    Ec = 4700*sqrt(fc′), #MPa
    fy = 420.0, #MPa
    Es = 210_000, #MPa
    ec_rebar = 5911.05e-9 #kgCO2e/mm3
    )

    d = input_vector[1] 
    bd_ratio = input_vector[2] 
    b = d*bd_ratio

    covering = 50.0 #mm 
    depth = -d-covering
    return depth
end


model = Model(beam_obj)
x0 = 1
addvar!(model,[0,0],[Inf,Inf],init = [1,1])
add_ineq_constraint!(model, beam_cons)
