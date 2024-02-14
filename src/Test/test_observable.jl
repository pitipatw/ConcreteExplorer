using Makie, GLMakie, CairoMakie
using DataFrames
GLMakie.activate!()

f1 = Figure(size = (500,500))
ax1 = Axis(f1[1,1])

x_all = rand(100)
y_all = rand(100)

df = DataFrame( [:x,:y] .=> [x_all, y_all])

s1 = IntervalSlider(f1[1,2], range= LinRange(0, 1, 1000), startvalues= (0,0), tellheight =false, horizontal = false)
s2 = IntervalSlider(f1[2,1], range= LinRange(0, 1, 1000), startvalues= (0,0), tellheight =true, horizontal = true)

#base plot 
scatter!(ax1, df[!, :x], df[!, :y], color = :grey)

points = rand(Point2f, 300)
colors = lift(s1.interval, s2.interval) do v1,v2
    map(points) do p
        (v2[1] < p[1] < v2[2]) && (v1[1] < p[2] < v1[2])
    end
end

scatter!(points, color = colors, colormap = [:gray90, :dodgerblue], strokewidth = 0)
