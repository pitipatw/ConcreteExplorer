using AsapSections
function find_Pu(c::ConcreteSection)
    return 0.65*0.85 * (c.fc′ * (c.geometry.area - sum(c.rebars.ast))) + sum(c.rebars.fy .* c.rebars.ast)
end
function find_Mu(c::ConcreteSection)
    area_req = sum(c.rebars.ast .* c.rebars.fy)/(c.fc′*0.85)
    a = depth_from_area(c.geometry,area_req, show_stats = false);
    β1 = clamp(0.85- 0.05*(c.fc′-28)/7, 0.65,0.85)
    c_ = a/β1
    d = c.geometry.ymax - c.geometry.ymin
    ϵs = 0.003 * (d - c_) / c_
    ϕ = clamp(0.65 + 0.25 * (ϵs - 0.002) / 0.003, 0.65, 0.90)
    mu = sum(c.rebars.ast .* c.rebars.fy .* (d - (a/2)))
    return  ϕ*mu
end

"""
Working on this...
"""
function find_fc′(Mu::Float64,section::AsapSections.PolygonalSection,rebars)
    #use bisection to find fc′.
    fc′_guess = 28
    step = 0.1
    tol = 1
    while tol > 1e-3
        c = ConcreteSection(fc′_guess, section, rebars)
        Mu_calc = find_Mu(c)
        tol = abs(Mu_calc-Mu)/Mu
        if count > 500
            println("Maximum Iteration Exceed, Please increase maxiter")
            return 0.0
        end
        fc′_guess = fc′_guess + step
    end
   return  fc′_guess
end


# function find_Vu(c::ConcreteSection)
#     println("Working in progress")
#     @assert 1 = 2
#     return 0.0
# end
