using Nonconvex, Zygote
# Nonconvex.@load MMA
Nonconvex.@load NLopt

include("beam_design.jl")

function beam_opt(fc′::Float64, ec::Float64)

    function beam_obj(input_vector::Vector{Float64};
        fc′::Float64 = fc′, ec::Float64 = ec)
        d = input_vector[1]
        bd_ratio = input_vector[2]
        return beam_design(d,bd_ratio,fc′,ec)[1]
    end

    """
    this function could be the same as obj,but with some other output 
    <Will work on this later>
    """
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
        h = d+covering
        area = b*h
        # p1 = [0.0 , 0.0]
        # p2 = [b, 0.0]
        # p3 = [b, h]
        # p4 = [0.0, h]

        # points = [p1, p2, p3, p4]
        # section = SolidSection(points)

        w_load = (1.2*wSDL + 1.6*wLL)*bay
        w_beamDL = 1.2*ρ_concrete*area
        w = w_load + w_beamDL
        Mdemand = w * span^2/8
        Vdemand = w * span/2

        C = 2*Mdemand/(0.85*fc′*b)
        B = -2*d
        
        con1 = -(B^2 - 4*C) #inside a root must be positive

        a = (-B - sqrt(B^2-4*C))/2
        a = a > 0 ? a : a = (-B + sqrt(B^2-4*C))/2

        con2 = a-d #value of a must be positive and less than d 

        as = Mdemand/(fy*(d-a/2))
        # c = a/0.85
        # ϵs = 0.003*(d-c)/c 
        # if ϵs < 0.005
        #     # println("NON ductile behavior")
        # end
        #0.0018 is from ACI (see the jupyter notebook)
        ρ_min = clamp(as/(area-as), 0.0018, Inf)

        #find ρ_max 
        c_max = 3/8*d
        a_max = 0.85*c_max
        as = Mdemand/(fy*(d-a_max/2))
        ρ_max = as/(area-as)

        #here, should we select the ρ_min always? 
        ρ_selected = (ρ_min+ρ_max)/2
        area_max = area*ρ_max
        area_min = area*ρ_min
        
        a_maximum = area_max*fy/(0.85*fc′*b)
        c_maximum = a_maximum/0.85
        ϵs_maximum = 0.003*(d-c_maximum)/c_maximum

        a_minimum = area_min*fy/(0.85*fc′*b)
        c_minimum = a_minimum/0.85
        ϵs_minimum = 0.003*(d-c_minimum)/c_minimum


        reinforcement_area = area*ρ_selected
        a = reinforcement_area*fy/(0.85*fc′*b)
        c = a/0.85
        ϵs = 0.003*(d-c)/c
        if ϵs < 0.005
            # println("Mean reinforcement is not ductile enough, using ρ_minimum instead")
            ρ_selected = ρ_min
            reinforcement_area = area*ρ_selected
            a = reinforcement_area*fy/(0.85*fc′*b)
            c = a/0.85
            ϵs = 0.003*(d-c)/c
        end

        # δ_max = 5*w*L^4/(384*E*I)

        #find E and I of the beam. 
        #we turn everything into a pure concrete with E = 4700sqrt(fc′) 
        as_transformed = Es/Ec*as 
        total_area = area-as+as_transformed
        #find a new neutral axis of the section 
        NA_concrete = -h/2# section.centroid[2]
        NA_rebar = -d
        NA_section = (area*NA_concrete+as_transformed*NA_rebar)/total_area

        @assert -d < NA_section < NA_concrete

        I_RCsection = 1/12*b*h^3 + area*(NA_section-NA_concrete)^2 + as_transformed*(NA_rebar-NA_section)^2

        # @assert section.Ix < I_RCsection

        δ_max = 5*w*span^4/(384*Ec*I_RCsection)
        #Deflection limit L/240 
        δ_limit =  span/240

        # @assert δ_max <= δ_

        #serviceability constraints 
        con3 = δ_max - δ_limit

        return con1, con2, con3 
    end 

    model = Model(beam_obj)
    lb = [200.0, 0.5]
    ub = [1000.0 ,1.0]
    
    #do a grid serach on the bouded box. 
    n = 5
    x1 = lb[1]:(ub[1]-lb[1]):ub[1]
    x2 = lb[1]:(ub[1]-lb[1]):ub[1]

    x0 = [1000.0, 0.75]

    println("Initial Guess", x0)
    println("Initial Obj: ", beam_obj(x0))


    addvar!(model, lb, ub, init = x0)
    add_ineq_constraint!(model, beam_cons)

    # alg = MMA02() # or MMA02()
    # options = MMAOptions(maxiter = 2000)
    alg = NLoptAlg(:LN_COBYLA)
    options = NLoptOptions()
    result = optimize(model, alg, x0, options = options, convcriteria = KKTCriteria() )
    println("*Opt design: " ,result.minimizer)
    println("*Opt obj: ",result.minimum)

    return result.minimizer, result.minimum
end