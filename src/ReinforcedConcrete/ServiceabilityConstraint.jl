module ServiceabilityConstraint
using Distributions
using Makie, GLMakie
#using PlotlyJS, DataFrames


#Assuming non crack, we will use I_g(gross Inertia) for calculation of beam defelction
# Setting up Equations
"""
Find the gross Inertia for a rectanglular section
Table 8-2 from EB070-07.pdf 
or Section 5, page 190 of RC Mechanics and Design 6th textbook
"""
function find_I_g(b, h)
    return (b * h^3) / 12
end

"""
Find the gross Inertia for a T section
"""
function find_I_g(b, b_w, h_f, h, y_t)
    return ((b - b_w) * h_f^3) / 12 + (b_w * h^3) / 12 + (b - b_w) * h_f * ((h - h_f) / (2 - y_t))^2 + (b_w * h * (y_t - (h / 2))^2)
end

"""
Find y_t(Distance from centroidal axis of gross section, neglecting reinforcement, to tension face)
Used in finding the gross Inertia for a T section
"""
function find_y_t(h, b, b_w, h_f)
    return h - (1 / 2) * ((b - b_w) * h_f^2 + (b_w * h^2)) / ((b - b_w) * h_f + (b_w * h))
end

"""
Modulus of elasticity for w_c(unitweight) between 1440 and 2560 kg/m^3
ACI 19.2.2.1
"""
function find_E_c(; w_c, fc′)
    return (w_c)^1.5 * 0.043 * sqrt(fc′) # in MPa
end

"""
Modulus of elasticity for normalweight concrete
ACI 19.2.2.1
"""
function find_E_c(fc′::Float64)
    return 4700.0 * sqrt(fc′) # in psi
end


"""
Ratio of A_s′(Area of compressive reinforcement) to bd
Calculate at midspan for simple and continuous spans and at the support for cantilevers 
ACI 24.2.4.1.2
"""
function find_ρ′(A_s′, b, d)
    return A_s′ / (b * d)
end

"""
Find time-dependent factor for sustained load
ACI Table 24.2.4.1.3
"""
function find_ξ(sustained_load_duration_in_months::Int64)
    if sustained_load_duration_in_months == 3
        return 1.0
    elseif sustained_load_duration_in_months == 6
        return 1.2
    elseif sustained_load_duration_in_months == 12
        return 1.4
    elseif sustained_load_duration_in_months >= 60
        return 2.0
    else
        return "Round to 3,6,12, 60 or more months"
    end
end

"""
Multiplier used for additional deflection due to long-term effects
"""
function find_λ_Δ(ξ, ρ′)
    return ξ / (1 + (50 * ρ′))
end

"""
Calculating Deflection of Simply Supported, Uniform Distributed Load
w in here is load
"""
function find_δ(w, x, L, E, I)
    return -(w * x) * (L^3 - (2 * L * x^2) + x^3) / (24 * E * I)
end

"""
Calculating the Maximum Deflection of simply supported, uniform Distribued load
which is at x= L/2
"""
function find_max_δ(w, L, E, I)
    return (5 * w * L^4) / (384 * E * I)
end

"""
Change ft to inches
"""
function convert_ft_to_in(ft)
    return ft * 12
end

### needed for stress???
"""
Find the modulus of rupture (f_r)
"""
function find_f_r(λ, fc′)
    return 0.62 * λ * sqrt(fc′)
end

"""
Obtain λ for lightweight
1.0 for normal weight
"""
function find_λ(w_c)
    if w_c <= 1600
        return 0.75
    elseif w_c > 1600 && w_c <= 2160
        return 0.0075 * w_c
    else
        return 1.0
    end
end
####
function main()

    # Define Initial values
    # Example 5-4 RC textbook
    Beam_Span_Length = 500 # in
    fc′ = rand(Uniform(2200, 4400), 30, 1) #4000.0 psi
    b_w = 12.0 #in
    h_f = rand(Uniform(0, 24), 30, 1) #in
    h = 24.0 #in
    d = h - 2.5 #in
    # d = Array{Float64,1}()
    # for h in h
    #     push!(d, h - 2.5) #in
    # end
    #σ = Fn / A
    M_u = 130000.0 #ft
    ϕ = 0.9 #no unit
    f_y = 60000.0 #psi
    j = 0.95 # assumption
    spacing_between_beams = 12.0 #ft
    dead_load = 2.5 #k/ft
    live_load = 1.5 #k/ft # kip/ square feet # kip= 1000 pound force
    f_s = (2 / 3) * f_y

    ## Deflection Limitation from Table 24.2.2 from ACI
    # Flat roofs (Not supporting/ not attached to elements likely to be damaged by large deflection)
    # l/180
    # Floor (Not supporting/ not attached to elements likely to be demaged by large deflection)
    # l/360
    # Roof or floors (supporting/ attached to nonstructural elements likely to be damaged by large deflection)
    # l/480
    # Roof or floors (supporting/ attached to nonstructural elements NOT likely to be damaged by large deflection)
    # l/240

    values_of_total_deflection_within_limits = Array{Float64,1}()
    values_of_total_deflection_out_of_limits = Array{Float64,1}()
    # beam_length_within_limits = Array{Float64,1}()
    # beam_length_out_of_limits = Array{Float64,1}()
    # dead_load_within_limits = Array{Float64,1}()
    # dead_load_out_of_limits = Array{Float64,1}()
    fc′_within_limits = Array{Float64,1}()
    fc′_out_of_limits = Array{Float64,1}()
    h_f_within_limits = Array{Float64,1}()
    h_f_out_of_limits = Array{Float64,1}()


    #Execution
    for each_fc′ in fc′
        for each_h_f in h_f
            # for each_length in Beam_Span_Length ##These two lines are for 3d plotting with varying beam length and deadloads
            #     for each_dead_load in D
            effective_width_b1 = convert_ft_to_in(Beam_Span_Length) / 4
            effective_width_b2 = b_w + (2 * 8 * each_h_f)
            effective_width_b3 = (convert_ft_to_in(spacing_between_beams))

            b = min(effective_width_b1, effective_width_b2, effective_width_b3)

            y_t = find_y_t(h, b, b_w, each_h_f)
            I_g = find_I_g(b, b_w, each_h_f, h, y_t)
            E = find_E_c(each_fc′)

            # Using basic combination 2 from 2.3.2 of ASCE standard
            # Using just dead loads first

            total_loads = (1.2 * dead_load * 1000 / 12) + (1.6 * live_load * 1000 / 12)  ## to change the units from kip to lbs and ft to in, multiple D and L by 1000/12
            #plotting_points[count, :] = [each_length, each_dead_load]
            max_δ = find_max_δ(total_loads, Beam_Span_Length, E, I_g)
            #Assuming that this is time dependent
            λ_Δ = find_λ_Δ(find_ξ(60), find_ρ′(0, b, d))
            #println("λ_Δ is ", λ_Δ)
            total_deflection = max_δ * λ_Δ
            #total_deflection = max_δ
            #println("total_deflection is ", total_deflection)

            deflection_limit = 25.4 * (Beam_Span_Length) / 180 ##bc each length is in in and this equation consider the length to be in mm
            #println("deflection limit is ", deflection_limit)

            if total_deflection < deflection_limit
                # push!(beam_length_within_limits, each_length)
                # push!(dead_load_within_limits, each_dead_load)
                push!(fc′_within_limits, each_fc′)
                push!(h_f_within_limits, each_h_f)
                push!(values_of_total_deflection_within_limits, total_deflection)
                #println(true)
            else
                # push!(beam_length_out_of_limits, each_length)
                # push!(dead_load_out_of_limits, each_dead_load)
                push!(fc′_out_of_limits, each_fc′)
                push!(h_f_out_of_limits, each_h_f)
                push!(values_of_total_deflection_out_of_limits, total_deflection)
                #println(false)
            end
        end
    end
    #visualization
    GLMakie.activate!()
    # tell julia to use GLMakie
    f = Figure(resolution=(1200, 800)) #initialize with resolution
    ax = Axis3(f[1, 1], xlabel="fc′[MPa] ", ylabel="h_f[in]", zlabel="max deflection")
    #initialize 3d axis with labels
    # scatter!(ax, beam_length_within_limits, dead_load_within_limits, values_of_total_deflection_within_limits, color=values_of_total_deflection_within_limits) #plot all data
    # scatter!(ax, beam_length_out_of_limits, dead_load_out_of_limits, values_of_total_deflection_out_of_limits, color="gray") #plot all data
    scatter!(ax, fc′_within_limits, h_f_within_limits, values_of_total_deflection_within_limits, color=values_of_total_deflection_within_limits) #plot all data
    scatter!(ax, fc′_out_of_limits, h_f_out_of_limits, values_of_total_deflection_out_of_limits, color="gray") #plot all data

    return f


    # Create a table of data
    # df = DataFrame(
    #     id=1:length(b_w_within_limits),
    #     b_w_values=b_w_within_limits,
    #     h_f_values=h_f_within_limits,
    #     h_values=h_within_limits,
    #     d_values=d_within_limits,
    #     deflection_values=values_of_total_deflection_within_limits
    # )

    # trace = parcoords(;
    #     line=attr(color=df.id),
    #     dimensions=[
    #         attr(range=[0, 40], label="b_w", values=df.b_w_values),
    #         attr(range=[0, 40], label="h_f", values=df.h_f_values),
    #         attr(range=[0, 40], label="h", values=df.h_values),
    #         attr(range=[0, 40], label="d", values=df.d_values),
    #         attr(range=[0, 500], label="deflection", values=df.deflection_values)
    #     ]
    # )
    # layout = Layout(
    #     title_text="Parallel Coordinates Plot",
    #     title_x=0.5,
    #     title_y=0,
    # )

    # df2 = DataFrame(
    #     id=1:length(b_w_out_of_limits),
    #     b_w_values=b_w_out_of_limits,
    #     h_f_values=h_f_out_of_limits,
    #     h_values=h_out_of_limits,
    #     d_values=d_out_of_limits,
    #     deflection_values=values_of_total_deflection_out_of_limits
    # )

    # trace2 = parcoords(;
    #     line=attr(color="red"),
    #     dimensions=[
    #         attr(range=[0, 40], label="b_w", values=df2.b_w_values),
    #         attr(range=[0, 40], label="h_f", values=df2.h_f_values),
    #         attr(range=[0, 40], label="h", values=df2.h_values),
    #         attr(range=[0, 40], label="d", values=df2.d_values),
    #         attr(range=[0, 500], label="deflection", values=df2.deflection_values)
    #     ]
    # )
    # layout2 = Layout(
    #     title_text="Parallel Coordinates Plot",
    #     title_x=0.5,
    #     title_y=0,
    # )


    # parallel_plot = plot(trace, layout)
    # parallel_plot2 = plot(trace2, layout2)
    # display(parallel_plot2)
    # display(parallel_plot)


end

end #module ServiceabilityConstraint

ServiceabilityConstraint.main()