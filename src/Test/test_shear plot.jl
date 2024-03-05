span::Float64 = 10_000.0
bay::Float64  = 10_000.0
wLL::Float64  = 2.40e-3#N/mm2
wSDL::Float64 = (24e-6*300)#2400kg/m3 * 300cm -> N/mm20
nV = 50


Vdemand = 1000* span/2
step_x = range(0,span,nV)
V(x) = abs(Vdemand-2*Vdemand/span*x)



using Makie, GLMakie
f1= Figure(size = (500,500))
a1 = Axis(f1[1,1])
scatter!(a1, step_x, V.(step_x))
