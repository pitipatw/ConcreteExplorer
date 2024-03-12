bd_ratio = 0.5
d = 600.0 #mm
b = bd_ratio*d #mm

p1 = [0.0 , 0.0]
p2 = [b, 0.0]
p3 = [b, -d]
p4 = [0.0, -d]

points = [p1, p2, p3, p4]
section = CompoundSection([SolidSection(points)])
section.solids
println(fieldnames(SolidSection))
f1 = Figure(size = (300,300)) 
a1 = Axis(f1[1,1], aspect = DataAspect() )
poly!(a1, section.solids[1].points)
f1


reinforcement_area = 4000
a = reinforcement_area/2 
r = sqrt(a/pi)

#create points on circle
covering = 50
offsety = -d+ covering 
offsetx = b/4
pts1 = [[r*cos(θ)+offsetx, r*sin(θ)+offsety] for θ in 0:2*pi/10:2pi]
pts2 = [[r*cos(θ)+b-offsetx, r*sin(θ)+offsety] for θ in 0:2*pi/10:2pi]

circle1 = SolidSection(pts1)
circle2 = SolidSection(pts2)
f1
poly!(a1, circle1.points, color = :red)
poly!(a1, circle2.points, color = :red)


rc_section1 = create_rc_section(section.solids[1], [circle1, circle2], d, 100.0)
rc_section1, serviceability = beam_design(46.0, 315.12e-9);
println(rc_section1.embodied_carbon)
f1 = draw(rc_section1)
