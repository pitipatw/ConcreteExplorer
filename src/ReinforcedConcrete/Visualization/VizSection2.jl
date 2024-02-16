using Makie, GLMakie
using kjlMakie
set_theme!(kjl_light)

function VizCatalog(catalog)
for s in 1:size(catalog)[1]
    # c is a concrete section 
    c = catalog[s, :Section]
    section_ID = catalog[s, :Section_ID]
    # rebar_section = getfield.(catalog[s,:Section],rebars) , 
    figure1 = Figure(resolution=(600, 600))
    xl = -25
    xu = 550
    yl = -550
    yu = 20
    ax1 = Axis(figure1[1, 1],# title="Section",
        aspect=DataAspect(),
        xticks=xl:50:xu, yticks=yl:50:yu,
        limits=(xl, xu, yl, yu))
    Label(figure1[1, 1, Top()], "Concrete section $section_ID\n Rebar section", valign = :bottom,
    font = :bold, fontsize = 30,
    padding = (0, 0, 5, 0))
    #plot section outline from AsapSections
    pts = c.geometry.points'
    #concate the first point to create a closed loop line
    pts = vcat(pts, pts[1, :]')
    lines!(ax1, pts[:, 1], pts[:, 2])

    #plot rebars from c.rebars
    bar_centers = hcat(c.rebars.x, c.rebars.y)

    # get circle points
    function circle_pts(r::Float64; n=50, base=[0.0, 0.0])
        return [r .* [cos(thet), sin(thet)] .+ base for thet in range(0, 2pi, n)]
    end

    #has to loop each rebar, but that's ok :) 
    for i = 1:length(c.rebars.d)
        #for visual only, we will sort the size first.
        center = bar_centers[i]
        d = c.rebars.d[i]
        x = c.rebars.x[i]
        y = c.rebars.y[i]
        r = d / 2
        pts = Point2f.(circle_pts(r, base=[c.rebars.x[i], c.rebars.y[i]]))
        lines!(ax1, pts, color=:blue, linewidth=2)
    end


    figure1
    save("src/ReinforcedConcrete/Sections/$s.png", figure1)
end

end

VizCatalog(catalog)
