"""
fc′ and fR1 and fR3 relationship

Available strengths

35/45

50/60

80/95

Available

"""

using CSV, DataFrames
println(@__DIR__)
fc_fiber  = CSV.read("src//fc_fiber.csv", DataFrame);

fc_fiber[!, :strength]

fc = collect(fc_fiber[!, :strength])
dosage = collect(fc_fiber[!, :dosage])