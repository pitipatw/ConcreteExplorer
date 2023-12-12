using Makie, GLMakie, kjlMakie

println("Sanity checking
Values from ASCE 7-10
table 4-1
Live load 
    for office use = 2.4 kN/m2
    for Computer use = 4.79 kN/m2
    Corridors = 4.79 kN/m2
****for residential
        Uninhabitable attics without storage 0.48
        Uninhabitable attics with storage 0.96
        Habitable attics and sleeping areas 1.44
        All other areas except stairs 1.92
        Private rooms and corridors serving them 1.92
        Public rooms and corridors serving them 4.79
    roof
        Garden 4.79
        ordinary flat, pitched, and curved roofs 0.96


in general: 1.92 kN/m2 should do the work.
Beware for 4.79 cases, since it's 2.5 times more than 1.92

Dead load 
    assume a 400x400 mm concrete section 
    Concrete density = 2400 kg/m3
    Dl = 0.4*0.4*2.4 = 0.384 kN/m - > 0.4 kN/m
")

#span goes from 2 meters to 12 meters 
#bay goes from 2 metrers to 8 meters

set_theme!(kjl_dark)
spans = Vector{Float64}()
bays = Vector{Float64}()
Mdemands_office = Vector{Float64}()
Mdemands_residential = Vector{Float64}()

set_spans = 2.0:0.5:12.0
set_bays  = 2:0.5:8

for span = set_spans
    for bay = set_bays
        w_office = 1.2 * 0.4 + 1.6 * (4.80 * bay)
        w_residential = 1.2 * 0.4 + 1.6 * (1.92 * bay)

        Mdemand_office = w_office * span^2 / 8
        Mdemand_residential = w_residential * span^2 / 8
        push!(spans, span)
        push!(bays, bay)
        push!(Mdemands_office, Mdemand_office)
        push!(Mdemands_residential, Mdemand_residential)
    end
end

f1 = Figure(resolution=(1000, 1000), title="Name")
ax1 = Axis(f1[1, 1], title="Beam spans vs Moment damands [Residential]",
    xlabel="Spans [m]", ylabel="Moment demands [kN.m]")
scatter!(ax1, spans, Mdemands_residential, markersize=10)

ax2 = Axis(f1[1, 2], title="Bays vs Moment damands [Residential] ",
    xlabel="Bays [m]", ylabel="Moment demands [kN.m]")
scatter!(ax2, bays, Mdemands_residential, markersize=10)

ax3 = Axis(f1[2, 1], title="Beam spans vs Moment damands [Office]",
    xlabel="Spans [m]", ylabel="Moment demands [kN.m]")
scatter!(ax3, spans, Mdemands_office, markersize=10, color=:green)

ax4 = Axis(f1[2, 2], title="Bays vs Moment damands [Office] ",
    xlabel="Bays [m]", ylabel="Moment demands [kN.m]")
scatter!(ax4, bays, Mdemands_office, markersize=10, color=:green)
f1

# save("span_bay_load.png", f1)