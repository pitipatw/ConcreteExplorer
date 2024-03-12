using LinearAlgebra
using AsapToolkit
using StaticArrays


"""
    makepixel(L::Real, t::Real, Lc::Real; n = 10)
<By Keith JL.>
L = length of pixel arm
t = thickness
Lc = straight region of pixel (length before arc)
n = number of discretizations for arc

Output : Vector of Vector of points [Vector{Vector{Float64}}]
"""
function make_pixel_geometry(L::Real, t::Real, Lc::Real; n=10)

    #constants
    θ = pi / 6
    ϕ = pi / 3
    psirange = range(0, ϕ, n)

    #origin
    p1 = [0.0, 0.0]

    # first set
    p2 = p1 .+ [0.0, -L]
    p2′ = p1 .+ L .* [cos(θ), sin(θ)]

    #second set
    p3 = p2 .+ [t, 0.0]
    p3′ = p2′ + t .* [cos(ϕ), -sin(ϕ)]

    #third set
    p4 = p3 .+ [0.0, Lc]
    p4′ = p3′ .+ Lc .* [-cos(θ), -sin(θ)]

    #arc
    v4 = p4′ .- p4

    #radius
    r = norm(v4) / cos(ϕ) / 2

    #arc center
    p5 = p4 .+ [r, 0.0]

    arcs = [p5 .+ r .* [-cos(ang), sin(ang)] for ang in psirange]

    points = [p1, p2, p3, arcs..., p3′, p2′]

    return points
end

#useful functions for rotating points
rotate_2d_about_origin(point::AbstractVector{<:Real}, angle::Float64) = [cos(angle) -sin(angle); sin(angle) cos(angle)] * point
rotate_2d_about_origin(point::Matrix{<:Real}, angle::Float64) = point * [cos(angle) -sin(angle); sin(angle) cos(angle)]
# rotate_2d_about_origin(point::Vector{Vector{Float64}}, angle::Float64) = point * [cos(angle) -sin(angle); sin(angle) cos(angle)] 
move_2d(point::Matrix{<:Real}, vector::Matrix{<:Real}) = point .+ vector

"""
By Keithjl
Make a Y (3-pieces) layup section
"""
function make_Y_layup_section(L::Real, t::Real, Lc::Real; n=10, offset=0.0)

    pixel = make_pixel_geometry(L, t, Lc; n=n)

    #offset from origin
    θ = pi / 6
    offset_vector = offset .* [cos(θ), -sin(θ)]

    #bottom right pixel
    right_pixel = [point + offset_vector for point in pixel]

    #top pixel
    top_pixel = rotate_2d_about_origin.(right_pixel, 2pi / 3)

    #bottom left pixel
    left_pixel = rotate_2d_about_origin.(top_pixel, 2pi / 3)

    sections = SolidSection.([right_pixel, top_pixel, left_pixel])

    return CompoundSection(sections)
end

"""
By Keithjl
Make a X2 (2-pieces X) layup section
"""
function make_X2_layup_section(L::Real, t::Real, Lc::Real; n=10, offset=0.0)
    pixel = make_pixel_geometry(L, t, Lc; n=n)

    #offset from origin
    θ = pi / 6

    offset_vector = offset .* [cos(θ), -sin(θ)]

    #base pixel
    base_pixel = [point + offset_vector for point in pixel]

    #right pixel 
    right_pixel = rotate_2d_about_origin.(base_pixel, pi / 6)

    #top pixel
    top_pixel = rotate_2d_about_origin.(right_pixel, pi / 2)

    #left pixel
    left_pixel = rotate_2d_about_origin.(top_pixel, pi / 2)

    #bottom pixel
    bottom_pixel = rotate_2d_about_origin.(left_pixel, pi / 2)

    distance = top_pixel[2][1] - right_pixel[end][1]

    right_pixel = [[point[1] + distance, point[2]] for point in right_pixel]
    top_pixel = [[point[1], point[2] + distance] for point in top_pixel]
    left_pixel = [[point[1] - distance, point[2]] for point in left_pixel]
    bottom_pixel = [[point[1], point[2] - distance] for point in bottom_pixel]

    sections = SolidSection.([top_pixel, bottom_pixel])
    return CompoundSection(sections)
    # return sections
end

"""
By Keithjl
Make a X4 (4-pieces X) layup section
"""
function make_X4_layup_section(L::Real, t::Real, Lc::Real; n=10, offset=0.0)

    pixel = make_pixel_geometry(L, t, Lc; n=n)

    #offset from origin
    θ = pi / 6
    offset_vector = offset .* [cos(θ), -sin(θ)]

    #base pixel
    base_pixel = [point + offset_vector for point in pixel]

    #right pixel 
    right_pixel = rotate_2d_about_origin.(base_pixel, pi / 6)

    #top pixel
    top_pixel = rotate_2d_about_origin.(right_pixel, pi / 2)

    #left pixel
    left_pixel = rotate_2d_about_origin.(top_pixel, pi / 2)

    #bottom pixel
    bottom_pixel = rotate_2d_about_origin.(left_pixel, pi / 2)

    distance = top_pixel[2][1] - right_pixel[end][1]

    right_pixel = [[point[1] + distance, point[2]] for point in right_pixel]
    top_pixel = [[point[1], point[2] + distance] for point in top_pixel]
    left_pixel = [[point[1] - distance, point[2]] for point in left_pixel]
    bottom_pixel = [[point[1], point[2] - distance] for point in bottom_pixel]

    sections = SolidSection.([right_pixel, top_pixel, left_pixel, bottom_pixel])
    return CompoundSection(sections)
    # return sections
end


viz = false
if viz == true
    include("../Visualization/plotpixel.jl")
    # viz_all_pixels()
end







#Will have to make a half pixel here.
#In the hope that I can just mirror that whole thing and make it faster? 
function halfpixel(L::Real, t::Real, Lc::Real; n=10)
    println("Hang in there. I'm working on it.")
end




# """
# Turn a vector of vector into matrix
# """
# function vecvec_to_matrix(vecvec)
#     dim1 = length(vecvec)
#     dim2 = length(vecvec[1])
#     my_array = zeros(Float64, dim1, dim2)
#     for i in 1:dim1
#         for j in 1:dim2
#             my_array[i, j] = vecvec[i][j]
#         end
#     end
#     return my_array
# end

# """
# fill a box that tighly confines pixel geometry with grid points, with a grid size dx, dy

# """
# function fillpoints(nodes::Matrix{Float64}, dx::Real, dy::Real)

#     pixel = make_pixel_geometry(L, t, Lc; n=n)

#     #offset from origin
#     θ = pi / 6
#     offset_vector = offset .* [cos(θ), -sin(θ)]

#     #base pixel
#     base_pixel = [point + offset_vector for point in pixel]

#     #right pixel 
#     right_pixel = rotate_2d_about_origin.(base_pixel, pi / 6)

#     #top pixel
#     top_pixel = rotate_2d_about_origin.(right_pixel, pi / 2)

#     #left pixel
#     left_pixel = rotate_2d_about_origin.(top_pixel, pi / 2)

#     #bottom pixel
#     bottom_pixel = rotate_2d_about_origin.(left_pixel, pi / 2)

#     distance = top_pixel[2][1] - right_pixel[end][1]

#     right_pixel = [[point[1] + distance, point[2]] for point in right_pixel]
#     top_pixel = [[point[1], point[2] + distance] for point in top_pixel]
#     left_pixel = [[point[1] - distance, point[2]] for point in left_pixel]
#     bottom_pixel = [[point[1], point[2] - distance] for point in bottom_pixel]

#     sections = SolidSection.([top_pixel, bottom_pixel])
#     # return CompoundSection(sections)
#     return sections
# end

