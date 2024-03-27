module Deflection

using AsapToolkit
using UnPack

#function to create pixelframe points
include("..//Geometry//pixelgeo.jl")

#pixelframe definition
include("pixelframe_element_property.jl")
include("pixelframe_element.jl")

#functions to calculate deflection
include("utilities.jl")
include("functions.jl")
include("interpolations.jl")

export Material,
    Properties

"""
function get_deflection(pixelframe_element::AbstractPixelFrameElement)
"""
function get_deflection(pixelframe_element::AbstractPixelFrameElement, applied_load::Real)
    #Unpack pixelframe_element into the root properties

    # @unpack fc′, Ec, Eps, fpy, em, es, em0, dps0, Ls, Ld, L, Aps, Atr, Itr, Zb, w, mg, fr, r, fpe, ϵpe, ϵce, Mdec, test = pixelframeelement

    # Mat = Material(fc′, Ec, Eps, fpy)
    # Sec = Section(em, es, em0, dps0, Ls, Ld, L, Aps, Atr, Itr, Zb)
    # f = Loads(w, mg, fr, r, fpe, ϵpe, ϵce, Mdec)
    @unpack material, properties = pixelframe_element

    return get_deflection(material, properties, applied_load)
end 

# function get_deflection(reinforcedconcrete_element::AbstractReinforcedConcreteElement)
#     @error "Function to find the deflection of a reinforced concrete is underway."
#     @unpack material, properties = reinforcedconcrete_element
#     return get_deflection_reinforcement_concrete
# end



"""
function get_deflection
Todo:
    add variable descriptions in the function.
    remove Mat, Sec, and f requirements.
    get_Deflection(pixelframeelement::PixelFrameElement, load::Float64)
Return mid-span deflection of a simply supported pixelframe element
"""
function get_deflection(material::Material,properties::Properties, applied_load::Real;
# function get_deflection(pixelframeelement::PixelFrameElement, load::Float64;
    loadstep::Float64=10.0)

    if loadstep > applied_load 
        loadstep = applied_load/100
    end

    @unpack fc′, Ec, Eps, fpy = material
    @unpack em, es, em0, dps0, Ls, Ld, L, Aps, Atr, Itr, Zb, w, moment_selfweight, fr, r, fpe, ϵpe, ϵce, moment_decompression = properties
    
    
    # we could do Mcr = 0 , becuase we crack at the begining anyway. 
    # but it happens to break the model instead. 
    # 0.00001 (very small) also breaks

    #Initializing parameters, so they exist globally.
    Ωc = 0
    c  = 0
    Ac_req = 0
    Lc = 0
    fc = 0.0
    δ_mid = 0
    δ_dev = 0
    Mi = 0
    # 4.448*8300.0 N
    P = 0.0:loadstep:load #P is the total load in the 2-point load test.
    M = P * Ls / 2.0 #given M inputs

    #Assumption
    Icr = Itr #Initial Icr. It will be calculated again later.
    Ie = Itr #Initial effective moment of inertia

    #Initial values assumption
    fps_assumed = fpe
    # dps = dps0

    #output containers and history for tracking variables.
    lenP    = length(P)
    fps_his = zeros(lenP)
    dps_his = zeros(lenP)
    Icr_his = zeros(lenP)
    Ie_his  = zeros(lenP)
    c_his   = zeros(lenP)
    dis_his = zeros(lenP)
    dis_dev_his = zeros(lenP)
    fc_his  = zeros(lenP)
    history = Vector{Vector{Float64}}([fps_his, dps_his, Icr_his, Ie_his, c_his, dis_his, dis_dev_his, fc_his])
    #Based on figure 7 on the paper. 

    Ω = getΩ(properties) 
    Mcr = getMcr(material, properties, Ω)

    for i in eachindex(M)
        Mi = M[i]
        if Mi <= Mcr 
            #DONE
            δ, fps, e = linearElasticUncrackRegime!(material, properties, Mi, fps_assumed, history, i)

        elseif Mcr < Mi < My 
            linearElasticCrackedRegime(Mi, pixelframe_element, )
        elseif My < Mi 
            ultimateRegime() 
        else 
            println("Invalid moment value")
        end 
        fps_assumed = fps
    end 

    return δ, fps

end





    for i in eachindex(M)
        Mi = M[i]
        Lc = clamp(L - 2 * Ls * Mcr / Mi, L - 2 * Ls, Inf)

        counter1 = 0
        conv1 = 1
        while conv1 > 1e-3
            counter1 += 1
            if counter1 > 1000
                println("Warning: 1st iteration did not converge")
                break
            end

            conv2 = 1
            counter2 = 0
            while conv2 > 1e-3
                # println("counter")
                counter2 += 1
                if counter2 > 1000
                    println("Warning: 2nd iteration did not converge")
                    break
                end
                #assume value of Itr and fps
                Ωc = getΩc(Ω, Icr, Lc, Sec)
                c = 10.0 #dummy
                conv_3 = 1
                counter_3 = 0
                while conv_3 > 1e-3
                    counter_3 += 1
                    if counter_3 > 1000
                        println("Warning: 3rd iteration did not converge")
                        break
                    end
                    #centroid of concrete area might not be at c/2
                    Ac_req = Mi / (dps - c / 2) / (0.85 * fc′)

                    # Ac_req = ps_force_i /0.85/fc′
                    if test
                        new_c = get_C(Ac_req, test=test)
                    else
                        new_c = get_C(pixelframeelement, Ac_req)
                    end

                    conv_3 = abs(new_c - c) / new_c
                    c = new_c
                end
                #calculate Icr
                if test
                    Icr_calc = get_Icrack(c, test=test)
                else
                    Icr_calc = get_Icrack(pixelframeelement, c)
                end

                conv2 = abs(Icr_calc - Icr) / Icr_calc

                Icr = Icr_calc

            end

            Ie = getIe(Mcr, Mdec, Mi, Icr, Itr)
            δ_mid, δ_dev, e = getDelta(Mat, Sec, f, Ie, Mi, em, fps)
            dps = dps0 - (δ_mid - δ_dev)
            fc = fps / Eps * c / (dps - c) + Mi / Itr * c  #concrete stress at top fiber
            fps_calc = getFps2(Mat, Sec, f, Ωc, c, dps, fc)
            conv1 = abs(fps_calc - fps) / fps
            fps = fps_calc
        end

        # δmid = getDeltamid()
        #record the history
        fps_history[i] = fps
        dps_history[i] = dps
        Icr_history[i] = Icr
        Ie_history[i] = Ie
        c_history[i] = c
        dis_history[i] = δ_mid
        fc_history[i] = fc
        dis_dev_history[i] = δ_dev
    end

    outputs = Dict("fps_history" => fps_history,
        "dps_history" => dps_history,
        "Icr_history" => Icr_history,
        "Ie_history" => Ie_history,
        "c_history" => c_history,
        "dis_history" => dis_history,
        "fc_history" => fc_history,
        "dis_dev_history" => dis_dev_history,
        "P" => P,
        "Mcr" => Mcr,
        "Mdec" => Mdec)
    if test
        return outputs
    else
        return dis_history, P
    end
end


end
