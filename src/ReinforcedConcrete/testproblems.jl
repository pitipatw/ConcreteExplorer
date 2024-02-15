# b = 300.0
# d = 350.0
# fc′ = 25.0
# fy = 500.0
# as = 1350.0
#unfactored mn = 210 kNm

b = 350
d = 450.0
fc′ = 32.0
fy = 500.0
as = 2480.0
#unfactored mn = 475 kNm
rebars = RebarSection([as], [fy], [b/2], [d], [0.])
p1 = [0. , 0.]
p2 = [b , 0.]
p3 = [b, -d]
p4 = [0., -d]
pts = [p1,p2,p3,p4]
section = SolidSection(pts)
c = ConcreteSection(fc′, section, rebars)
M = find_Mu(c)/1e6

