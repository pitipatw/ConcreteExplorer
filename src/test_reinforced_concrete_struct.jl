using AsapSections

a = collect([ 0. 0. ; 1. 0. ; 1. 1. ; 0. 1.]')
b = a .+ [1. 0.]'

function circle_pts(r::Float64; n=50, base=[0.0, 0.0])
    return [r .* [cos(thet), sin(thet)] .+ base for thet in range(0, 2pi, n)]
end

section = CompoundSection(SolidSection.([a,b]))

offset = [-100 0]
circle_pts.([10.0 20.0], base = offset)
section.solids[1].points


concrete_section_pts = collect([0.0 0.0 ; 300.0 0 ; 300.0 -400.0 ; 0.0 -400]')
rebar_radius = [10.0 30.0 10.0]
rebar_pos = [100.0 150.0 200.0 ; -350.0 -300.0 -350.0]
concrete_carbon = 280.0

create_rc_section(concrete_section_pts, rebar_radius, rebar_pos, concrete_carbon)