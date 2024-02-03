using UnPack

include("definition.jl") #Will have to merge this with the normal stuff.
include("functions.jl")
include("interpolations.jl")
include("structuralelement.jl")

"""
Todo:
    add variable descriptions in the function.
    remove Mat, Sec, and f requirements.
    get_Deflection(pixelframeelement::PixelFrameElement, maximum_load::Float64)
Return mid-span deflection of a simply supported pixelframe element
"""
function get_Deflection(pixelframeelement::PixelFrameElement, maximum_load::Float64;
    loadstep::Float64=10.0)::Dict{String,Any}

    @unpack fc′,Ec,Eps,fpy,em,es,em0,dps0,Ls,Ld,L,Aps,Atr,Itr,Zb,w,mg,fr,r,fpe,ϵpe,ϵce,Mdec,test = pixelframeelement
    @show test

    Mat = Material(fc′, Ec, Eps, fpy)
    Sec = Section(em, es, em0, dps0, Ls, Ld, L, Aps, Atr, Itr, Zb)
    f   = Loads(w, mg, fr, r, fpe, ϵpe, ϵce, Mdec)

    #Could do
    # if test
    #     get_Deflection(pixelframeelement)
    # else
    #     get_Deflection(pixelframeelement)
    # end

    #Initializing parameters
    Ωc = 0
    c  = 0
    Ac_req  = 0 
    Lc = 0
    fc = 0.0
    δ_mid = 0
    δ_dev = 0
    Mi = 0
    # 4.448*8300.0 N
    P = 0.0:loadstep:maximum_load
    M = P*Ls/2.0 #given M inputs

    #Assumption
    Icr = Itr #Initial Icr. It will be calculated again later.
    Ie = Itr #Initial effective moment of inertia

    #Initial values
    fps = fpe 
    dps = dps0 
    Ω =  getOmega(Sec)
    #we could do Mcr = 0 , becuase we crack at the begining anyway. 
    #but it happens to break the model instead. 0.00001 (very small) also breaks
    @show Mcr = getMcr(Mat, Sec, f, Ω)


    #output containers
    lenP = length(P)
    fps_history     = zeros(lenP)
    dps_history     = zeros(lenP)
    Icr_history     = zeros(lenP)
    Ie_history      = zeros(lenP)
    c_history       = zeros(lenP)
    dis_history     = zeros(lenP)
    dis_dev_history = zeros(lenP)
    fc_history      = zeros(lenP)

    #Based on figure 7 on the paper.

    conv1 = 1
    counter1 = 0
    counter2 = 0
    for i in eachindex(M)
        Mi = M[i] 
        Lc = getLc(Sec,Mcr,Mi)

        counter1 = 0
        conv1 = 1
        while conv1 > 1e-3
            counter1 += 1 
            if counter1 > 1000
                println("Warning: 1st iteration did not converge")
                break
            end
            #assume value of Itr and fps

            conv2 = 1
            counter2 = 0
            while conv2 > 1e-3
                # println("counter")
                counter2 += 1 
                if counter2 > 1000
                    println("Warning: 2nd iteration did not converge")
                    break
                end
                Ωc = getΩc(Ω, Icr, Lc, Sec)
                # ps_force_i = Aps*fps

                c = 10.0 #dummy
                conv_c = 1 
                counter_c = 0 
                while conv_c > 1e-3 
                    counter_c += 1
                    if counter_c > 1000
                        println("Warning: 3rd iteration did not converge")
                        break
                    end
                    #centroid of concrete area might not be at c/2
                    Ac_req = Mi/(dps-c/2)/(0.85*fc′)
                
                    # Ac_req = ps_force_i /0.85/fc′
                    new_c = get_C(Ac_req, test = test)
                    conv_c = abs(new_c - c)/new_c
                    c = new_c
                end
                #calculate Icr
                Icr_calc = get_Icrack(c, test = test)

                conv2 = abs(Icr_calc - Icr)/Icr_calc

                Icr = Icr_calc
                
            end
        
            Ie = getIe(Mcr, Mdec, Mi, Icr, Itr)
            δ_mid, δ_dev , e  = getDelta(Mat, Sec, f, Ie, Mi, em,fps)
            dps = dps0 - (δ_mid - δ_dev)
            fc = fps/Eps*c/(dps-c) + Mi/Itr*c  #concrete stress at top fiber
            fps_calc = getFps2(Mat, Sec, f , Ωc, c, dps, fc)
            conv1 = abs(fps_calc - fps) / fps
            fps = fps_calc
        end

        # δmid = getDeltamid()
        #record the history
        fps_history[i] = fps
        dps_history[i] = dps
        Icr_history[i] = Icr
        Ie_history[i]  = Ie
        c_history[i]   = c
        dis_history[i] = δ_mid
        fc_history[i]  = fc
        dis_dev_history[i] = δ_dev
    end

    outputs = Dict("fps_history" => fps_history,
    "dps_history" => dps_history,
    "Icr_history" => Icr_history,
    "Ie_history" => Ie_history,
    "c_history" =>  c_history,
    "dis_history" => dis_history,
    "fc_history" =>  fc_history,
    "dis_dev_history" => dis_dev_history,
    "P" => P,
    "Mcr" => Mcr, 
    "Mdec" => Mdec)
    
    return outputs
end 

