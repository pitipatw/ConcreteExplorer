include("reinforced_concrete.jl")

concrete_section_pts = collect([0.0 0.0 ; 300.0 0 ; 300.0 -400.0 ; 0.0 -400]')
rebar_radius = [15.0 ,15.0, 15.0]
rebar_pos = [100.0 150.0 200.0 ; -350.0 -330.0 -350.0]
concrete_carbon = 280.0

[circle_pts(rebar_radius[i], base = rebar_pos[:,i]) for i in eachindex(rebar_radius)]
rebar_sections = CompoundSection(SolidSection.([circle_pts(rebar_radius[i], base = rebar_pos[:,i]) for i in eachindex(rebar_radius)]))

s1 = create_rc_section(concrete_section_pts, rebar_radius, rebar_pos, concrete_carbon)
draw(s1)