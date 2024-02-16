module ParallelComputing
using PlotlyJS, DataFrames

# Initializing the domain range
β_1 = 0.85
b = Array{Float64}(120.0:100:1200.0) # mm
h = Array{Float64}(120.0:100:1270.0) # mm
fc′ = Array{Float64}(28.0:2:50.0) # MPa
As = Array{Float64}(100.0:50:500.0) # sq mm
f_y = 413 #MPa
d = Array{Float64,1}()
for each_h in h
    push!(d, each_h - 63.5) # mm
end


# Setting up a empty matrix for all the combinations
total_combination_length = length(b) * length(fc′) * length(As) * length(d)
total_combination_matrix = Matrix{Float64}(undef, total_combination_length, 5)
count = 0

# Filling in the matrix with the combinations 
for each_b in b
    for each_fc′ in fc′
        for each_As in As
            for each_d in d
                global count += 1
                total_combination_matrix[count, :] = [each_b, each_fc′, each_As, each_d, f_y]
            end

        end
    end
end

## Writing functions
"""
Equation to determine the depth of the compression stress block, a
This equation can be used to find A_s in 5-16
"""
function find_a(A_s, f_y, fc′, b)
    return (A_s * f_y) / (0.85 * fc′ * b)
end

"""
This is taken straight from Eq.5-15 * ϕ in RC textbook (Example 5-4)
"""
function calculate_flexural_strength(b, fc′, As, d, f_y)
    a = find_a(As, f_y, fc′, b)
    c = a / β_1
    if c < (3 / 8) * d
        ϕ = 0.9
    else
        ϕ = 0.65
    end
    flexual_strength = ϕ * As * f_y * (d - a / 2)
    return flexual_strength
end

## this is where all the execution will happen
function main(num_of_iteration)

    ## algorithm to randomly pick a (number of) row(s) in the matrix created above
    ## and calculate the flexural strength and save them

    #create a random number generator and also set up empty Arrays and Matrix to save values
    random_number_generator = rand(1:total_combination_length, num_of_iteration, 1)
    results = Array{Float64,1}()
    values_for_plotting = Matrix(undef, num_of_iteration, 5)

    # distributing the work of calculating flexural strength of chosen rows to multiple threads
    # and saving them
    Threads.@threads for i = 1:num_of_iteration
        values_for_plotting[i, :] = total_combination_matrix[random_number_generator[i], :]
        push!(results, calculate_flexural_strength(values_for_plotting[i, 1], values_for_plotting[i, 2], values_for_plotting[i, 3], values_for_plotting[i, 4], values_for_plotting[i, 5]))
    end

    #Create a table of data
    df = DataFrame(
        id=1:num_of_iteration,
        b_values=values_for_plotting[:, 1],
        fc′_values=values_for_plotting[:, 2],
        As_values=values_for_plotting[:, 3],
        d_values=values_for_plotting[:, 4],
        f_y_values=values_for_plotting[:, 5],
        flexural_strength_values=results
    )

    trace = parcoords(;
        line=attr(color=df.flexural_strength_values),
        dimensions=[
            attr(range=[10, 1300], label="b[mm]", values=df.b_values),
            attr(range=[0, 100], label="fc′[MPa]", values=df.fc′_values),
            attr(range=[50, 550], label="As[sq mm]", values=df.As_values),
            attr(range=[10, 1300], label="d[mm]", values=df.d_values),
            attr(range=[10, 1300], label="f_y[MPa]", values=df.f_y_values),
            attr(range=[1e4, 2.3e8], label="flexural strength", values=df.flexural_strength_values)
        ]
    )
    layout = Layout(
        title_text="Parallel Coordinates Plot",
        title_x=0.5,
        title_y=0,
    )

    ##Create a parallel plot
    parallel_plot = plot(trace, layout)
    return parallel_plot
end
end


## Put in the number of iteration wanted
ParallelComputing.main(70)

