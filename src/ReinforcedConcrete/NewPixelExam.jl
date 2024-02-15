# module NewPixelExam
# end
# # #define packages that will be used
# # using Makie, GLMakie
# # using StatsBase


# # #setting up equations
# # """
# # Change ft to inches
# # """
# # function convert_ft_to_in(ft)
# #     return ft * 12
# # end

# # """
# # Equation for determining the required steel area in a singly reinforced section
# # 5-16 RC textbook
# # """
# # function find_A_s(M_u::Float64, ϕ::Float64, f_y::Float64, d::Float64, a::Float64)
# #     return M_u / (ϕ * f_y * (d - (a / 2)))
# # end
# # """
# # 5-16 Rc textbook equation but with different variables
# # """
# # function find_A_s(; M_u::Float64, ϕ::Float64, f_y::Float64, j::Float64, d::Float64)
# #     return M_u / (ϕ * f_y * (j * d))
# # end
# # """
# # Equation to determine the depth of the compression stress block, a
# # This equation can be used to find A_s in 5-16
# # """
# # function find_a(A_s, f_y, fc′, b)
# #     return (A_s * f_y) / (0.85 * fc′ * b)
# # end

# # """
# # Equation for minimum A_s
# # """
# # function find_A_smin(fc′::Float64, b_w::Float64, d::Float64, f_y::Float64)
# #     return (3 * sqrt(fc′) * b_w * d) / f_y
# # end

# # ### In this example, we will be assuming it's a single layer reinforcement ###
# # ### Meaning, it is not preferred to do combinations of diff sizes bar###
# # ### so we will keep at 2 bars of same size combinations
# # function main()

# #     # initial known values

# #     M_u = 130000.0 #ft
# #     ϕ = 0.9 #no unit
# #     f_y = 60000.0 #psi
# #     j = 0.95 # assumption
# #     height = Array{Float64}(10.0:0.5:50.0)
# #     h = sample(height, 60) #in #Height or the overall thickness of member
# #     d = Array{Float64,1}()
# #     for h in h
# #         push!(d, h - 2.5)
# #     end
# #     println(d)
# #     b_w = 12.0 #in #web width, diameter of circular section
# #     beam_span_length = 28.0 #ft
# #     h_f = 6.0 #in
# #     spacing_between_beams = 12.0 #ft
# #     fc′ = 4000.0 #psi


# #     #creating a dictionary for two same size bars combinations'areas
# #     two_bars_area = Dict("No.3" => 0.22, "No.4" => 0.4, "No.5" => 0.62,
# #         "No.6" => 0.88, "No.7" => 1.2, "No.8" => 1.58, "No.9" => 2,
# #         "No.10" => 2.54, "No.11" => 3.12, "No.14" => 4.5, "No.18" => 8)

# #     #creating a dictionary for bar size and their radius
# #     bars_radius = Dict("No.3" => 0.1875, "No.4" => 0.25, "No.5" => 0.3125,
# #         "No.6" => 0.375, "No.7" => 0.4375, "No.8" => 0.5, "No.9" => 0.564,
# #         "No.10" => 0.635, "No.11" => 0.705, "No.14" => 0.8465, "No.18" => 1.1285)


# #     #creating an array to store points to plot


# #     #execution
# #     values_of_d = Array{Float64,1}()
# #     values_of_actual_A_s = Array{Float64,1}()
# #     strength = Array{Float64,1}()
# #     for d in d
# #         push!(values_of_d, d)
# #         A_s = find_A_s(M_u=convert_ft_to_in(M_u),
# #             ϕ=ϕ, f_y=f_y, j=j, d=d)

# #         A_smin = find_A_smin(fc′, b_w, d, f_y)

# #         println("Initial A_s =", A_s)
# #         println("A_s minimum = ", A_smin)

# #         # first checkin point
# #         if A_smin < A_s
# #             # Section 4-8 of RC textbook, ACI 8.12.2
# #             # limits for the effective width of the compression flange
# #             effective_width_b1 = convert_ft_to_in(beam_span_length) / 4
# #             effective_width_b2 = b_w + (2 * 8 * h_f)
# #             effective_width_b3 = (convert_ft_to_in(spacing_between_beams))

# #             b = min(effective_width_b1, effective_width_b2, effective_width_b3)

# #             println(b)

# #             #find a using the smallest b
# #             a = find_a(A_s, f_y, fc′, b)
# #             println("a= ", a)

# #             # find the better A_s value to find the bars size combination
# #             better_A_s_value = find_A_s(convert_ft_to_in(M_u), ϕ, f_y, d, a)
# #             println(better_A_s_value)
# #             push!(values_of_actual_A_s, better_A_s_value)
# #             #interpolation to find the best possible bars combo
# #             #or should I do this manually???


# #             #Calculating the spacing limit once the number of bar and size are decided

# #             # bars selection
# #             number_of_bars = 2
# #             f_s = 40000 #psi
# #             β1 = 0.85 #this is for the concrete compressive strength fc′ of 4000 psi
# #             c_c = 2 #in
# #             A_s_final = 0.0
# #             bar_size = "No.0"
# #             #this new A_s is due to the combination of bar selected, in this case, 2 No.8 bars
# #             smallest_diff_btwn_required_As_and_combination_As = 99999
# #             for (BarSize, Area) in two_bars_area
# #                 current_diff = Area - better_A_s_value
# #                 if current_diff < smallest_diff_btwn_required_As_and_combination_As && current_diff > 0
# #                     smallest_diff_btwn_required_As_and_combination_As = current_diff
# #                     A_s_final = Area
# #                     bar_size = BarSize
# #                 end
# #             end
# #             println("A_s_final = ", A_s_final)
# #             println("bar_size = ", bar_size)
# #             # side to bar center distance depends on bar size
# #             side_to_bar_center_dist = 2.065 + get(bars_radius, bar_size, 0)

# #             # write the center to center spacing limit function
# #             s1 = 15 * (40000 / f_s) - (2.5 * c_c)
# #             s2 = 12 * (40000 / f_s)

# #             s = min(s1, s2)
# #             println(s)

# #             actual_s = b_w - (number_of_bars * side_to_bar_center_dist)
# #             println("actual s = ", actual_s)

# #             if actual_s < s
# #                 # if the code passes both check points, proceed to this required strength check
# #                 a_final = find_a(A_s_final, f_y, fc′, b)

# #                 c = a_final / β1
# #                 println("c = ", c)

# #                 # checkpoint 3
# #                 # check if c is less than 3/8 of d, to ensure this is a tension-controlled section

# #                 if c < (3 * d) / 8
# #                     # check again to see if jd(d-a/2) is larger than asuumed value
# #                     # which is 0.95d in this case and how much percent larger
# #                     final_jd = d - (a_final / 2)
# #                     initial_jd = 0.95 * d
# #                     percent_larger_btwn_final_and_initial_jd = 100 * (final_jd - initial_jd) / abs(initial_jd)

# #                     println("Final jd is larger than initial assumption by ", percent_larger_btwn_final_and_initial_jd, "%")

# #                     #final final check
# #                     M_n = A_s_final * f_y * (final_jd)
# #                     strength_of_final_section_design = ϕ * M_n
# #                     if strength_of_final_section_design >= M_u
# #                         push!(strength, strength_of_final_section_design)
# #                     else
# #                         push!(strength, 0.0)
# #                     end
# #                 else
# #                     println(false)
# #                 end


# #             else
# #                 println(false)
# #             end
# #         else
# #             push!(values_of_actual_A_s, 0.0)
# #             push!(strength, 0.0)
# #         end
# #     end

# #     #visualization
# #     GLMakie.activate!()
# #     # tell julia to use GLMakie
# #     f = Figure(resolution=(1200, 800)) #initialize with resolution
# #     ax = Axis3(f[1, 1], xlabel="Effective Depth [In]", ylabel="Area of steel required[in^2]", zlabel="Strength of final section design [k-ft]")
# #     #initialize 3d axis with labels
# #     scatter!(ax, values_of_d, values_of_actual_A_s, strength, color=strength) #plot all data

# #     return f

# # end

# # end # module NewPixelExam

# # NewPixelExam.main()