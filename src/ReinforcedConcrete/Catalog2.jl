using AsapSections
using DataFrames
include("Definitions.jl")
include("inversefcprime.jl")

"""
    Design equation: 
        Mᵤ= ϕAₛfy(d-a/2)
     0 
    We enforce tension-controlled section
    i.e., ϵc = 0.003 and ϵs = 0.005 -> ϕ = 0.9
    We will have 
        c = 3/8d
    a = β1c , where β1 is a function of fc′
    
    Mᵤ= ϕAₛfy(d-a/2) becomes 
      = 0.9Aₛfy( d - β1(3/8)d/2 )

    Then, we can impose a sc

"""
function find_Mu2(c)
    d = c.rebars.y
    c_ =  3/8*d  
    β1 = clamp(0.85- 0.05*(c.fc′-28)/7, 0.65,0.85)
    a = β1*c_
    return 0.9*sum(c.rebars.ast .* c.rebars.fy .* (d .- (a/2)))
end


# function get_catalog2()
#The design space
set_Mu = 1e6:1e6:1e7
set_ρs = 0.01:0.005:0.08 #reinforcement ratio 
set_d = 200.:25.:500.
set_bd_ratio = 0.5:0.05:1 #b/d ratio b = bd_ratio*d
fy = 420.
covering = 40. #ACI318M-19 Table 20.5.1.3.1, Not exposed to weather or in contact with ground

catalog = Dict()
Cs = Vector{ConcreteSection}()
Ps = Vector{Float64}()
Ms = Vector{Float64}()
GWPs = Vector{Float64}()
Section_IDs = Vector{Int64}()

# section_ID = 0
count = 0 

for Mu in set_Mu
for ρs in set_ρs
    for d in set_d
        for bd_ratio in set_bd_ratio
            # section_ID = section_ID + 1
            count = count +1 
            #Create a rectangular section
            b = bd_ratio*d
            w = b
            h = d + covering
            p1 = [0. , 0.]
            p2 = [w , 0.]
            p3 = [w, -h]
            p4 = [0., -h]
            pts = [p1,p2,p3,p4]
            section = SolidSection(pts)
            as = ρs*b*d
            #NO Spacing check
            #Only 1 rebar (see as a whole)
            areas = [as]
            fys = [fy]
            xs = [w/2]
            ys = [d]
            ds = [sqrt(as/pi)]
            rebars = RebarSection(areas, fys, xs, ys, ds)
            
            ff = define_function(as, fy, b, d, Mu)
            try 
                fc′ = find_fc′(ff) 
            catch
                fc′ = 0
            end

            if fc′ > 55 || fc′< 28
                fc′ = 0.
                continue
            end

            # if ϵs > 0.005 #only work on this case.
            c = ConcreteSection(fc′, section, rebars)
            #Calculate the capacity.
            # P = find_Pu(c)
            M = find_Mu2(c)
            # V = find_Vn(c)
            push!(Cs,c)
            # push!(Ps,P)
            push!(Ms,M)
            push!(GWPs, c.gwp)
            push!(Section_IDs, section_ID)
        end
    end
end
end

    # after got all of the catalog, compared them 
catalog = DataFrame(id = 1:count, Section = Cs, Mu= Ms, Gwp = GWPs, Section_ID = Section_IDs)
println("Done catalog")
#     return catalog
# end

catalog = get_catalog2();
