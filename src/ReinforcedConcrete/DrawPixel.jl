module DrawPixel
using LinearAlgebra
using Makie, GLMakie
using GeometryBasics

"""
makepixel(L::Real, t::Real, Lc::Real; n = 10)
L = length of pixel arm
t = thickness
Lc = straight region of pixel (length before arc)
n = number of discretizations for arc
"""
function makepixel(L::Real, t::Real, Lc::Real; n=100)
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
    return points, p5, r
end

function fullpixel(L::Real, t::Real, Lc::Real; n=100)
    g1 = makepixel(L, t, Lc, n=n)
    ptx1 = [i[1] for i in g1[1]]
    pty1 = [i[2] for i in g1[1]]
    #remove first point (0.0)
    ptx1 = ptx1[2:end]
    pty1 = pty1[2:end]
    # ptx = vcat(ptx1, -ptx1)
    # pty = vcat(pty1, pty1)

    ptx = ptx1
    pty = pty1
    nodes = [ptx pty]
    #rotate to the top
    newpoints1 = Matrix{Float64}(undef, size(nodes)[1], 2)
    # draw a full pixelframe section
    for i = 1:size(nodes)[1]
        x = nodes[i, 1]
        y = nodes[i, 2]
        r = sqrt(x^2 + y^2)
        θ = atand(y / x)
        newθ = θ + 120.0
        newx = r * cosd(newθ)
        newy = r * sind(newθ)
        newpoints1[i, :] = [newx, newy]
    end

    #rotate to the side (flip)
    newpoints2 = Matrix{Float64}(undef, size(nodes)[1], 2)
    # draw a full pixelframe section
    for i = 1:size(nodes)[1]
        x = nodes[i, 1]
        y = nodes[i, 2]
        r = sqrt(x^2 + y^2)
        θ = atand(y / x)

        newθ = θ + 240.0
        newx = r * cosd(newθ)
        newy = r * sind(newθ)
        newpoints2[i, :] = [newx, newy]
    end
    newnodes = vcat(nodes, newpoints1, newpoints2)
    return newnodes
end



function halfpixel(L::Real, t::Real, Lc::Real; n=100)
    g1 = makepixel(L, t, Lc, n=n)
    ptx1 = [i[1] for i in g1[1]]
    pty1 = [i[2] for i in g1[1]]
    #remove first point (0.0)
    ptx1 = ptx1[2:end]
    pty1 = pty1[2:end]
    # ptx = vcat(ptx1, -ptx1)
    # pty = vcat(pty1, pty1)

    ptx = ptx1
    pty = pty1
    nodes = [ptx pty]
    #rotate to the top
    newpoints1 = Matrix{Float64}(undef, size(nodes)[1], 2)
    # draw a full pixelframe section
    for i = 1:size(nodes)[1]
        x = nodes[i, 1]
        y = nodes[i, 2]
        r = sqrt(x^2 + y^2)
        θ = atand(y / x)
        newθ = θ + 120.0
        newx = r * cosd(newθ)
        if newx > 0 #we only need the right half to mirror
            newy = r * sind(newθ)
            newpoints1[i, :] = [newx, newy]
        else
            newpoints1[i, :] = [0, 0]
        end
    end


    newnodes = vcat(nodes, newpoints1,)
    return newnodes
end

function main()
    points = halfpixel(150, 20, 10)
    #println(points)
    negative_points = -1 * (points)
    new_points = Matrix{Float64}(undef, size(points)[1], 2)
    #println(new_points)
    for i in 1:size(points)[1]
        x = negative_points[i, 1]
        y = points[i, 2]
        new_points[i, :] = [x, y]
    end

    f1 = Figure(resolution=(800, 800))
    ax1 = Axis(f1[1, 1], xlabel="x", ylabel="y", aspect=DataAspect())#, aspect = DataAspect(), xgrid = false, ygrid = false)
    scatter!(ax1, points[:, 1], points[:, 2], color=:green)
    scatter!(ax1, new_points[:, 1], new_points[:, 2], color=:green)
    lines(new_points[1, :], new_points[2, :])

    return f1
end

# function main()
#     points = makepixel(150, 10, 10)
#     f1 = Figure(resolution=(600, 600))
#     ax1 = Axis(f1[1, 1], aspect=DataAspect(), xlabel="x", ylabel="y")
#     poly!(Point2.(points))
#     return f1
# end

end

DrawPixel.main()