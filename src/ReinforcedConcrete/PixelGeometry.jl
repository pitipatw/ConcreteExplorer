module PixelGeometry
using LinearAlgebra
using GeometryBasics
using Makie, GLMakie

"""
    makepixel(L::Real, t::Real, Lc::Real; n = 10)
L = length of pixel arm
t = thickness
Lc = straight region of pixel (length before arc)
n = number of discretizations for arc
"""
function makepixel(L::Real, t::Real, Lc::Real; n=10)
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
    points = [p1, p2, p3, p4, arcs..., p4′, p3′, p2′]
    return points
end

"""
rotate a set of points by desired degrees 
"""
function rotate_the_points(points::Vector, Degrees::Float64)
    # change vector(vector(Float64)) into a matrix of float 64 
    # to perform matrix operation
    point_matrix = Matrix{Float64}(undef, size(points)[1], 2)
    for row in 1:size(point_matrix)[1]
        point_matrix[row, :] = points[row]
    end
    rotation_matrix = [cosd(Degrees) -sind(Degrees); sind(Degrees) cosd(Degrees)]
    NewPoints = Array{Array{Float64,1},1}()
    # rotate each point in the matrix and store them as new points 
    for row in 1:size(point_matrix)[1]
        rotated_points = rotation_matrix * point_matrix[row, :]
        x = rotated_points[1]
        y = rotated_points[2]
        push!(NewPoints, [x, y])
    end
    return NewPoints
end

"""
Output a full pixel with desired orientation degree
"""
function execute(rotation_degree::Float64)
    extracted_points = makepixel(150, 10, 10)
    # create a full pixel with desired degree of orientation
    points = rotate_the_points(extracted_points, rotation_degree)
    new_points = rotate_the_points(extracted_points, rotation_degree + 120.0)
    new_points1 = rotate_the_points(extracted_points, rotation_degree + 240.0)
    # plot the full pixel
    f1 = Figure(resolution=(600, 600))
    ax1 = Axis(f1[1, 1], aspect=DataAspect(), xlabel="x", ylabel="y")
    poly!(Point2.(points))
    poly!(Point2.(new_points))
    poly!(Point2.(new_points1))
    return f1
end

end

PixelGeometry.execute(45.0)