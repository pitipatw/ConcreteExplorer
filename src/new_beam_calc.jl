#New RC calculation for optimization scheme

function section_calc(d::Float64, bd_ratio::Float64, fc′::Float64, ec_concrete::Float64, as::Float64, asv::Float64;
    fy::Float64 = 420.0,
    Es::Float64 = 20_000.0, 
    covering::Float64 = 50.0,
    )

    b = d*bd_ratio
    area = b*d
    h = d+covering
    a = as*fy/(0.85*fc′*b)
    β1 = clamp(0.85 - 0.05 * (fc′ - 28.0) / 7, 0.65, 0.85)
    c = a/β1
    ϵs = (d-c)/c*0.003
    ϕ = clamp(0.65 + 0.25 * (ϵ - 0.002) / 0.003, 0.65, 0.90)

    Mcapacity = ϕ*0.85*fc′*b*d*(d-a/2)

    Ec = 4700*sqrt(fc′)
    as_transformed = Es/Ec*as 
    total_area = area-as+as_transformed
    #find a new neutral axis of the section 
    NA_concrete = -h/2 #section.centroid[2]
    NA_rebar = -d
    NA_section = (area*NA_concrete+as_transformed*NA_rebar)/total_area


    I_RCsection = 1/12*b*h^3 + area*(NA_section-NA_concrete)^2 + as_transformed*(NA_rebar-NA_section)^2
    
    δ_limit =  span/240


    #shear design 

    vc = 0.17*sqrt(fc′)
    vs = s*sldkf*sdlfkj



end 


function constrains() 
    ϵs, Mcapacity, δ_max = section_calc() 
    con1 = ϵs - 0.005 
    con2 = Mcapacity - Mdemand 
    con3 = δ_limit - δ_max
    return con1, con2, con3
end
 