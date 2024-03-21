using AsapToolkit
using Makie, GLMakie
#PM diagram 
include("..//Functions//Geometry//pixelgeo.jl")

L = 205; t = 35; Lc = 35 
section = make_Y_layup_section(L,t,Lc)

fc′ = 28
fy = 420
w = 200
d = 400
as = 100

ϵc = 0.003

ϵs = ϵc*(d-c)/c

#pure compression 
pn = 0.85*fc′*(w*d-as) + fy*as
pu = 0.65*0.8*pn

#zero tension at the tendon.
c = d
β1 = clamp(0.85 - 0.05 * (fc′ - 28.0) / 7, 0.65, 0.85)

a = β1*c

pn = 0.85*fc′*w*a 
mn = 0.85*fc′*w*a*(d/2-a/2) 

pu = 0.65*pn
mn = 0.65*mn
#
#fs = 0.5fy

#fs = fy (ϵs = 0.002)

#ϵs = 0.005 (max)

#pure bending.