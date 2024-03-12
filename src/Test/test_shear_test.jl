#check moment capacity of the shear test. 
include("Functions//Geometry//pixelgeo.jl")
include("Functions//capacities.jl")



fc′ = 35.0
Δx = 20.0 #mm gap between loading point and support.
post_tension_force = 142 #kN
tendon_area =  314.0 #mm2
stress_on_tendon = post_tension_force/tendon_area
max_axial = post_tension_force #kN
max_shear_interface = 0.3*post_tension_force #kN
applied_force = 2*max_shear_interface

max_moment = applied_force*Δx/1000 /2  #kNm
max_shear_applied = applied_force/2 #kN #should have been this number.


L  = 205.0
Lc = 30.0
t  = 35.0
section = make_Y_layup_section(L,t,Lc)
section_σ = force/section.area
v_max = get_Vu(section, fc′, 2.0, 2.0, tendon_area, stress_on_tendon, L)

get_Vu(section, fc′, 2.0, 2.0, tendon_area, 0.0, L)

#axial capacity 
axial_capacity = 0.85*section.area*fc′/1000 #kN

#moment capacity
# mm3 * N/mm2 -> Nmm, Nmm/10^6 = kNm
moment_capacity = section.Sx*fc′/10e6
section.Sx
section.Ix/section.ymax
section.Ix/section.ymin
#shear capacity 
shear_capacity = 0.17*sqrt(fc′)*section.area/1000

axial_ratio = max_axial/axial_capacity
moment_ratio = max_moment/moment_capacity 
shear_ratio = max_shear/shear_capacity