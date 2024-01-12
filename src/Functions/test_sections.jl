pts1 = copy(transpose([0. 0.; 10. 0.; 10. 10.; 0. 10.].*10))
pts2 = copy(transpose([[0. 10.]; [10. 10.]; [5. 30.]].*10))

s1 = SolidSection(pts1)
s2 = SolidSection(pts2)

c1 = CompoundSection([s1,s2])

f1 = Figure(resolution = (100,100))
ax1= Axis(f1[1,1])

lines!(ax1, c1.solids[1].points)
lines!(ax1, c1.solids[2].points)


s1.area
s2.area
c1.area