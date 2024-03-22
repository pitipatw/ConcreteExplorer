section = make_Y_layup_section(205,35,30)
ct = section.ymax - section.centroid[2] 
cb = section.centroid[2] - section.ymin

st = section.Ix/ct 
sb = section.Ix/cb


L = 5000
I = section.Ix
area = section.area
pi*16*2
fc′ = 40
Ec = 4700*sqrt(fc′)
ft =  0.25*sqrt(fc′)
fc = 0.6*fc′
w = 2400*9.81*area/1e9 #N/mm

δmax =  5/384*w*L^4/(Ec*I)

δlimit = 5000/240

#Axial post tensioned to prevent the 5m beam, under selfweight from open up.
P = w*L^2*area/(8*sb) #N
#draped profile, counter the deflection from the selfweight.
P2 = 6*5*w*L^4*24/(384*1500*(15000*5000-4*1500^2))
P3 = δmax*24*Ec*I*6/(0.3*(3-0.36)*5000^3) #N 
P4 = w*L^2/(8*sb*(1/area + 350/(sb))) #N 


P5 = 3e7/(sb*(1/area+205/sb)) #N
δmax_2nd =  5/384*w*5500^4/(Ec*I)
P6 = δmax_2nd*8*Ec*I/(5500^2 *205)