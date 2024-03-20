include("Functions/gen_demand.jl")

demands = gen_demand()

open("src//20_03_test_demands.json", "w") do f
    JSON.print(f, demands)
end

#Generate demands for results in the paper

#8 meters beams with 4 meters total bay width.
#live load 
#       1.92 kN/m2, 
#dead load
#       beam weight 2400 kg/m3 (concrete section)
#       Approx section size = 0.1 m3 ->  make_Y_layup_section(300,35,35).area/1e6 ~ 0.1 m3
#       10cm thick floor 240 (2400*0.1) kg/m2

bay_width = 4 #m
LL = 1.92*bay_width #kN/m

assumed_section_sizes = 0.1
DL = 2.4*9.81*assumed_section_sizes + 2.4*9.81*0.1*bay_width # kN/m

#beam length of 8m.

L = 8000 #mm 

#Primary

