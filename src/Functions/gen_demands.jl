#create beams demands from the varying 4-point load test.
using DataFrames
using JSON
P = 10:10:1000 #kN?
Types = ["primary", "secondary", "columns"]
span_length = 6000.0
S_idx = Int.(1:span_length/500)
E_idx = 1:length(P);
fields = [:e_idx, :s_idx, :pu, :mu, :vu, :ec_max, :type]

function PMV(p,x,t,L1,L2)
    @show total_L = 2*L1+L2
    if t == "columns"
            Pout = p
            Vout = 0 
            Mout = 0 
    elseif (t=="primary") || (t=="secondary")
            Pout = 0
            if x < L1 
                Vout = p/2 
                Mout = x*p/2
            elseif L1<= x <= L1+L2
                Mout = L1*p/2
                Vout = 0
            elseif L1+L2 < x <= total_L
                Vout = p/2
                Mout = L1*p/2 - (x-L1-L2)*p/2
            else
                println("x Error")
            end
    else
        println("Error, Invalid section type")
    end

    return Pout, Vout, Mout 
end



function get_demands(p::Real, t::String;
    L1::Real=2000, L2::Real=2000,
    step::Real=500)
    
    total_L = 2*L1+L2
    n_sections = Int(total_L/step)
    demand_P = zeros(n_sections)
    demand_V = zeros(n_sections)
    demand_M = zeros(n_sections)

    #We will look at the middle cut of each segment.
    midsegpos = 250:500:total_L
    for ix in eachindex(midsegpos)
        x = midsegpos[ix]
        lx = x-249
        ux = x+249
        lip, liv, lim = PMV(p,lx, t, L1,L2)
        uip, uiv, uim = PMV(p,ux, t, L1,L2)
        
        demand_P[ix] = maximum([lip, uip])
        demand_V[ix] = maximum([liv, uip])
        demand_M[ix] = maximum([lim, uip])
    end

    return demand_P, demand_V, demand_M 
end

# a = get_demands(100,"primary")

demands = DataFrame([name => [] for name in fields])
for ip in eachindex(P)
    p = P[ip]
    for t in Types
        Pdemand,Vdemand,Mdemand = get_demands(p, t)
        #each P, V, and M is a full demand on the beam.
        #We have to distribute them onto each index of the output json (DataFrame).
        for s_idx in S_idx 
            # fields = [:e_idx, :s_idx, :pu, :mu, :vu, :ec_max, :type]
            entry = [ip,s_idx, Pdemand[s_idx], Mdemand[s_idx], Vdemand[s_idx], 0, t]
            # @show entry = Dict(fields .=> [ip,s_idx, P[s_idx], M[s_idx], V[s_idx], Inf, t])
            push!(demands,entry)
        end
    end
end

demands[!, :idx] = 1:size(demands)[1]

open("src//test_demands.json", "w") do f
    JSON.print(f, demands)
end


loaded_demands = JSON.parsefile(joinpath("src//test_demands.json"));
