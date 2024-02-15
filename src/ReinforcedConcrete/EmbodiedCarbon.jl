module EmbodiedCarbon
## The script to calculate embodied carbon, forces capacity, and serviceability constrain

#define packages that will be used
using Distributions
using Makie, GLMakie
using GeometryBasics
#using StatsPlots

#setting up equations
"""
Change ft to inches
"""
function convert_ft_to_in(ft)
    return ft * 12
end

"""
Equation for determining the required steel area in a singly reinforced section
5-16 RC textbook
"""
function find_A_s(M_u::Float64, ϕ::Float64, f_y::Float64, d::Float64, a::Float64)
    return M_u / (ϕ * f_y * (d - (a / 2)))
end
"""
5-16 Rc textbook equation but with different variables
"""
function find_A_s(; M_u::Float64, ϕ::Float64, f_y::Float64, j::Float64, d::Float64)
    return M_u / (ϕ * f_y * (j * d))
end

"""
Equation for minimum A_s
"""
function find_A_smin(fc′::Float64, b_w::Float64, d::Float64, f_y::Float64)
    return (3 * sqrt(fc′) * b_w * d) / f_y
end

function find_A_smin(c::ConcreteSection)
    fc′ = c.fc′
    w = c.w
    d = c.h-50.0 #height - Covering
    fy = c.rebars.fy[1] #let's get the first one for now. Usually we use the same fy
    return (3 * sqrt(fc′) * w * d) / fy
end

"""
Equation to determine the depth of the compression stress block, a
This equation can be used to find A_s in 5-16
"""
function find_a(A_s, f_y, fc′, b)
    return (A_s * f_y) / (0.85 * fc′ * b)
end

function find_a(c::ConcreteSection)
    as = sum(c.rebars.ast)
    fy = c.rebars.fy[1]
    fc′ = c.fc′
    area = as*fy/(0.85*fc′)
    section = c.geometry
    return depth_from_area(section, area)
end


### Clean this up!!!

"""
input the calculated A_s value and match it with discrete rebar combinations
return the rebar combination and the steel area
"""
function give_minimum_A_s_rebar_combination(A_s::Float64)

    # create a dictionary with index to bar area to allow easy execution of for loop
    index_to_bar_area = Dict(
                            1 => 0.11, 
                            2 => 0.2,
                            3 => 0.31, 
                            4 => 0.44, 
                            5 => 0.6,
                            6 => 0.79, 
                            7 => 1.0, 
                            8 => 1.27,
                            9 => 1.56,
                            10 => 2.25,
                            11 => 4,
                            )

    # match the index to actual bar numbers
    index_to_bar_num = Dict(
                            1 => "No.3",
                            2 => "No.4",
                            3 => "No.5",
                            4 => "No.6",
                            5 => "No.7",
                            6 => "No.8",
                            7 => "No.9",
                            8 => "No.10",
                            9 => "No.11",
                            10 => "No.14",
                            11 => "No.18",
                            )

    # first make a dictionary of using 1 bar only 
    bar_combination_and_area = Dict(
                            "No.3" => 0.11, 
                            "No.4" => 0.2, 
                            "No.5" => 0.31,
                            "No.6" => 0.44, 
                            "No.7" => 0.6,
                            "No.8" => 0.79,
                            "No.9" => 1,
                            "No.10" => 1.27,
                            "No.11" => 1.56, 
                            "No.14" => 2.25, 
                            "No.18" => 4,
                            )

    @assert length(index_to_bar_num) == length(index_to_bar_area)
    @assert length(index_to_bar_num) == length(bar_combination_and_area)

    # add the two same bars combination into the dictionary above
    # doesn't need to calculate all combinations because the reinforcement needs to be symmetry
    # meaning if there's only 2 bars, they have to be the same size
    for each_rebar in 1:11
        total_steel_area = index_to_bar_area[each_rebar] * 2
        push!(bar_combination_and_area, string("2 ", index_to_bar_num[each_rebar]) => total_steel_area)
    end

    # add three bars combination into the dictionary 
    # only need to account for combinations with two same digits such as 122, since it needs symmetry 
    for first_rebar in 1:11
        for second_rebar in 1:11
            total_steel_area = index_to_bar_area[first_rebar] + (2 * index_to_bar_area[second_rebar])
            push!(bar_combination_and_area, string(index_to_bar_num[first_rebar], " & 2 ", index_to_bar_num[second_rebar]) => total_steel_area)
        end
    end

    
    #precalc.
    A_s_final = 0.0
    bar_size = "No.0"
    smallest_diff_btwn_required_A_s_and_combination_A_s = 99999

    for (BarSize, ComboArea) in bar_combination_and_area
        current_diff = ComboArea - A_s
        if current_diff < smallest_diff_btwn_required_A_s_and_combination_A_s && current_diff > 0
            smallest_diff_btwn_required_A_s_and_combination_A_s = current_diff
            A_s_final = ComboArea
            bar_size = BarSize
        end
    end
    return [bar_size, A_s_final]
end

"""
input the calculated A_s value and give all possible rebar combination
return the rebar combination and the steel area
"""
function give_all_possible_A_s_rebar_combination(A_s::Float64)
##repak
    # create a dictionary with index to bar area to allow easy execution of for loop
    index_to_bar_area = Dict(
                            1 => 0.11,
                            2 => 0.2,
                            3 => 0.31,
                            4 => 0.44, 
                            5 => 0.6,
                            6 => 0.79,
                            7 => 1,
                            8 => 1.27,
                            9 => 1.56,
                            10 => 2.25,
                            11 => 4
                            )

    # match the index to actual bar numbers
    index_to_bar_num = Dict(
                            1 => "No.3",
                            2 => "No.4",
                            3 => "No.5",
                            4 => "No.6",
                            5 => "No.7",
                            6 => "No.8",
                            7 => "No.9",
                            8 => "No.10",
                            9 => "No.11", 
                            10 => "No.14",
                            11 => "No.18"
                            )

    # first make a dictionary of using 1 bar only 
    bar_combination_and_area = Dict(
                                "No.3" => 0.11,
                                "No.4" => 0.2,
                                "No.5" => 0.31,
                                "No.6" => 0.44,
                                "No.7" => 0.6,
                                "No.8" => 0.79,
                                "No.9" => 1, 
                                "No.10" => 1.27,
                                "No.11" => 1.56,
                                "No.14" => 2.25, 
                                "No.18" => 4
                                )

    # add the two same bars combination into the dictionary above
    # doesn't need to calculate all combinations because the reinforcement needs to be symmetry
    # meaning if there's only 2 bars, they have to be the same size
    for each_rebar in 1:11
        total_steel_area = index_to_bar_area[each_rebar] * 2
        push!(bar_combination_and_area, string("2 ", index_to_bar_num[each_rebar]) => total_steel_area)
    end

    # add three bars combination into the dictionary 
    # only need to account for combinations with two same digits such as 122, since it needs symmetry 
    for first_rebar in 1:11
        for second_rebar in 1:11
            total_steel_area = index_to_bar_area[first_rebar] + (2 * index_to_bar_area[second_rebar])
            push!(bar_combination_and_area, string(index_to_bar_num[first_rebar], " & 2 ", index_to_bar_num[second_rebar]) => total_steel_area)
        end
    end

    all_A_s_combo_bigger_than_input_A_s_list = Array{Float64,1}()
    all_bar_size_combo_bigger_than_input_A_s_list = Array{String,1}()

    for (BarSize, ComboArea) in bar_combination_and_area
        each_pair_difference = ComboArea - A_s
        if each_pair_difference >= 0
            push!(all_A_s_combo_bigger_than_input_A_s_list, ComboArea)
            push!(all_bar_size_combo_bigger_than_input_A_s_list, BarSize)
        end
    end
    # println(all_bar_size_combo_bigger_than_input_A_s_list)
    return all_A_s_combo_bigger_than_input_A_s_list
    #return all_A_s_combo_bigger_than_input_A_s_list, all_bar_size_combo_bigger_than_input_A_s_list
end


"""
P0 22.4.2.2
Nominal Axial Compressive Strength
fc′= specified compressive strength of concrete
Ast= Area of reinforcement/ steel (A_s)
Ag= Area of concrete cross section
fy= specified yield strength for nonprestressed reinforcement
"""
function find_P0(fc′::Number, Ag::Float64, Ast::Float64, fy::Float64)
    return 0.85 * (fc′ * (Ag - Ast)) + (fy * Ast)
end

function findPn(c::ConcreteSection)
    return 0.85 * (c.fc′ * (c.ag - sum(c.rebars.ast))) + sum(c.rebar.fy .* c.rebar.ast)
end

"""
Nominal Moment Capacity
"""
function find_Mn(A_s, f_y, d, a)
    return A_s * f_y * (d - (a / 2))
end

function find_Mn(c::ConcreteSection)
    #calculate a 
    println("Calc a (compression depth) please")
    return sum(c.rebar.as .* c.rebar.fy .* (c.rebar.d .- (a/2)))
end

"""
Shear Capacity
"""
function find_shear_capacity()
    return phi * V_n
end

function find_Vn(c)
    #assume stirrup spacing_between_beam

    return 0.0

end



"""
Calculate the Embodied Carbon per unit
For concrete: ~2400 kg/cu.m
For rebar: ~7850 kg/cu.m
"""
function CalculateEmbodiedCarbon(b, h, length)
    MaterialVolume = b * h * length
    EmbodiedCarbon = coefficient * MaterialVolume
    return EmbodiedCarbon
end

## meeting update
# geometry will be fixed, input the steel position and steel size range
# this equation is derived from example 5-4 in RC textbook
function find_fc′(A_s::Float64, f_y::Float64, d::Float64, b::Float64, j::Float64)
    fc′ = (A_s * f_y) / (-2 * 0.85 * d * b * (j - 1))
    return fc′
end

"""
Calculating the Maximum Deflection of simply supported, uniform Distribued load
which is at x= L/2
"""
function find_max_δ(w, L, E, I)
    return (5 * w * L^4) / (384 * E * I)
end

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

function main()
    ## Define initial value

    # Later, this function will execute these inputs
    # 1.Section width (b)
    # 2.Section height (h)
    #   a.Reinforcement position (d) -> could be fixed to h - 50mm (so we don’t have too much variables)
    # 3.Number of reinforcement bars (let’s do 1 to 8 bars, at most)
    # 4.Reinforcement bars area (let’s assume same area for all bars, using the discrete bar sizes you already have in your script)
    # 5.Beam length (l)(This is crucial for serviceability analysis and weight comparison)

    # Define All Initial Values
    b = 5.0 #in
    b_w = 12.0 #in
    h_f = 12.0 #in
    h = 5.0 #in    
    d = h - 0.05 #in
    j = rand(Uniform(0, 5), 20, 1) # no unit
    f_y = 60000.0 #psi
    ϕ = 0.9 #no unit
    l = 10 # beam length in inch
    M_u = rand(Uniform(100000, 200000), 10, 1) #ft
    spacing_between_beams = 12.0 #ft
    dead_load = 2.5 #k/ft
    live_load = 1.5 #k/ft # kip/ square feet # kip= 1000 pound force

    # Creating the storage place for the variables
    j_within_limits = Array{Float64,1}()
    j_out_of_limits = Array{Float64,1}()
    M_u_within_limits = Array{Float64,1}()
    M_u_out_of_limits = Array{Float64,1}()
    A_s_within_limits = Array{Float64,1}()
    A_s_out_of_limits = Array{Float64,1}()
    fc′_within_limits = Array{Float64,1}()
    fc′_out_of_limits = Array{Float64,1}()
    values_of_total_deflection_within_limits = Array{Float64,1}()
    values_of_total_deflection_out_of_limits = Array{Float64,1}()
    j_when_fc′_is_negative = Array{Float64,1}()
    M_u_when_fc′_is_negative = Array{Float64,1}()
    A_s_when_fc′_is_negative = Array{Float64,1}()
    negative_fc′ = Array{Float64,1}()

    # With varying j and M_u, get multiples A_s and for each A_s, compare it with the combination of bar we have
    # and output all the A_s greater than or equal to the input A_s
    # lastly, find fc′ for each A_s and store all the variables
    for each_j in j
        for each_M_u in M_u
            A_s = find_A_s(M_u=convert_ft_to_in(each_M_u),
                ϕ=ϕ, f_y=f_y, j=each_j, d=d)
            available_A_s = give_all_possible_A_s_rebar_combination(A_s)
            for each_A_s in available_A_s
                fc′ = find_fc′(each_A_s, f_y, d, b, each_j)

                if fc′ > 0
                    ## Serviceability
                    y_t = find_y_t(h, b, b_w, h_f)
                    I_g = find_I_g(b, b_w, h_f, h, y_t)
                    E = find_E_c(fc′)

                    # Using basic combination 2 from 2.3.2 of ASCE standard
                    # Using just dead loads first

                    total_loads = (1.2 * dead_load * 1000 / 12) + (1.6 * live_load * 1000 / 12)  ## to change the units from kip to lbs and ft to in, multiple D and L by 1000/12
                    max_δ = find_max_δ(total_loads, l, E, I_g)
                    #Assuming that this is time dependent
                    λ_Δ = find_λ_Δ(find_ξ(60), find_ρ′(0, b, d))

                    total_deflection = max_δ * λ_Δ
                    #total_deflection = max_δ

                    deflection_limit = 25.4 * (l) / 180 ##bc each length is in in and this equation consider the length to be in mm
                    if total_deflection < deflection_limit
                        push!(fc′_within_limits, fc′)
                        push!(A_s_within_limits, each_A_s)
                        push!(j_within_limits, each_j)
                        push!(M_u_within_limits, each_M_u)
                        push!(values_of_total_deflection_within_limits, total_deflection)
                    else
                        println(false)
                        push!(fc′_out_of_limits, fc′)
                        push!(A_s_out_of_limits, each_A_s)
                        push!(j_out_of_limits, each_j)
                        push!(M_u_out_of_limits, each_M_u)
                        push!(values_of_total_deflection_out_of_limits, total_deflection)
                    end
                elseif fc′ < 0
                    push!(j_when_fc′_is_negative, each_j)
                    push!(M_u_when_fc′_is_negative, each_M_u)
                    push!(A_s_when_fc′_is_negative, each_A_s)
                    push!(negative_fc′, fc′)
                end
            end
        end
    end


    #visualization
    GLMakie.activate!()
    # tell julia to use GLMakie

    f1 = Figure(resolution=(1200, 800)) #initialize with resolution
    ax1 = Axis3(f1[1, 1], xlabel="j", ylabel="A_s", zlabel="Compressive Strength of Concrete")
    #initialize 3d axis with labels
    scatter!(ax1, j_within_limits, A_s_within_limits, fc′_within_limits, color=fc′_within_limits) #plot all data
    scatter!(ax1, j_when_fc′_is_negative, A_s_when_fc′_is_negative, negative_fc′, color="gray") #plot all data

    return f1

end
end


EmbodiedCarbon.main()