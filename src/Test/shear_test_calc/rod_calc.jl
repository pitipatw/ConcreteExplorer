include("pixelgeo.jl")
println("The Pixelframe we are testing has L,t,Lc of 205 35 30 mm")
println("This bring the area to:")
L = 205; t =35; Lc = 30;
section = make_Y_layup_section(L,t,Lc)
@show area = section.area
println("($area mm2)")

println("For post tensioning forces")
println("The tendon dia is 15.2mm -> 140 mm2 crosssection area. ")

println("We Load maximum at at 70% increment of 1860 MPa")
@show max_tendon_σ = 0.7*1860
@show max_f = max_tendon_σ*140*2
@show concrete_σ = max_f/area

println("Max stress $max_tendon_σ MPa")
println("maximum_force $max_f N")
println("Stress on concrete $concrete_σ MPa")

println("The one that we have size 18537.69 mm2")
old_area = 18537.69
test_force = concrete_σ * old_area
test_force_kN = test_force/1000

area/old_area




