include("Functions/gen_demand.jl")

demands = gen_demand()

open("src//16_03_test_demands.json", "w") do f
    JSON.print(f, demands)
end