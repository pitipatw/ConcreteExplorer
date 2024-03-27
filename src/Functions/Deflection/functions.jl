"""
Figure 6 Linear elastic uncracked regime.
***for the first guess of fps, fpe is preferred.
This function (workflow) includes.
(2) Tendon stress from the linear elastic uncracked section analysis.

if the tendon profile and the loaded points arent follow fig 4,\n 
equation (7) through (14) needs to be applied.

"""
function linearElasticUncrackRegime!(material::Material, properties::Properties, Mi::Float64, fps_assumed::Float64, history::Vector{Vector{Float64}},i::Int64)
    # @unpack Mat, Sec, F, type = element
    @unpack fc′, Ec, Eps, fpy, fr = material
    @unpack Aps, Atr, Itr, Zb, em, es, Ls, Ld, L, w, Mg, fr, r, fpe = properties
    @unpack em, es,θ, em0, dps0, L, Ls, Ld, Aps, Atr,Itr,r,Zb,w,moment_selfweight,fpe,moment_decompression,ϵpe, ϵce = properties

    #(6)
    Ω = getΩ(properties)

    if Ld < Ls 
        K1 = Ls/L - 1
        K2 = Ls/Ls(Ld/L)^2 - (Ls/L)^2 
    elseif Ld ≥ Ls
        K1 = Ld/L - 1
        K2 = 0
    end

    tol = 1 
    counter = 0
    while tol > 1e-3
        counter += 1 
        if counter > 1000
            println("getfps exceeds its limit")
            return  0
        end

        if type == "primary"
            top_e = em +Mi*L^2/(6*Ec*Itr)*(3*Ld/L*(-K1) - 3/4 - K2)
            bot_e = 1 - fps_assumed*Aps/(Ec*Itr)*(L^2/8 - L*Ld/2 + Ld^2/2)
            e =  top_e/bot_e 
        elseif type == "secondary" || type == "columns"
            k1 = 1/72 ; k2 = 1/72 ; 
            Δ = k1*Mi*L^2/(Ec*Itr) - k2*fps_assumed*Aps*e*L^2/(Ec*Itr)
            e = em-Δ
        else 
            println("Invalid section type")
            return 0.0
        end

        new_fps = fpe + Ω*Mi*e/(Itr*Ec/Eps + Aps*(r^2+e^2)*Ω)
        tol = abs(new_fps-fps_assumed)/fps_assumed
        fps_assumed = new_fps
    end

    fps = fps_assumed

    if type == "primary" 
        δ = Mi * L^2 / (6 * Ec * Itr) * (3 / 4 - (Ls / L)^2) - fps * Aps / (Ec * Itr) * (e * L^2 / 8 - (e - es) * Ls^2 / 6)
    elseif type == "secondary" || type == "columns"
        k3 = 1/8; 
        δ = 23*Mi*L^2/(216*Ec*Itr) - k3*fps*Aps*e*L^2/(Ec*Itr)
    end
    # history :[fps_his, dps_his, Icr_his, Ie_his, c_his, dis_his, dis_dev_his, fc_his]

    history[1][i] = fps 
    history[2][i] = e
    history[6][i] = δ

    return δ, fps, e
end

"""
Figure 7 Linear elastic uncracked regime.
***for the first guess of fps, fpe is preferred.
This function (workflow) includes.
(2) Tendon stress from the linear elastic uncracked section analysis.

if the tendon profile and the loaded points arent follow fig 4,\n 
equation (7) through (14) needs to be applied.

"""
function linearElasticCrackedRegime(M::Float64, fps_assumed::Float64, Icr_assumed::Float64, element::Element)
    @unpack Mat, Sec, F,geometry, type = element
    @unpack fc′, Ec, Eps, fpy = Mat
    @unpack Aps, Atr, Itr, Zb, em, es, Ls, Ld, L = Sec
    @unpack w, Mg, fr, r, fpe, Mdec = F
    
    Ω =  getΩ(Sec) #(6)
    Mcr = getMcr(Mat, Sec, F, Ω) 
    Lc = clamp(L-2*Ls*Mcr/M, L-2*Ls, Inf) #(20)
    
    #establish counter for the outermost loop.
    counter1 = 0
    conv1 = 1
    fps = fps_assumed #define fps here so it exists at this scope.
    while conv1 > 1e-3
        counter1 += 1 
        if counter1 > 1000
            println("Warning: 1st iteration did not converge")
            break
        end
            
        conv2 = 1
        counter2 = 0

        Icr = Icr_assumed
        c_assumed = 10.0 #assumed c.
        c = c_assumed #define here so it exist at this scope.
        while conv2 > 1e-3
            # println("counter")
            counter2 += 1 
            if counter2 > 1000
                println("Warning: 2nd iteration did not converge")
                break
            end

            Ωc = getΩc(Ω, Icr_assumed, Lc, Sec) # (21)
        
            #now, we try to calculate the compression depth c.
            # c = getC(c_assumed,Ωc,element)
            #calculate Icrack based on c. 

            if test #this is based on the interpolation csv file.
                Icr = get_Icrack(c, test = test)
            else 
                Icr = get_Icrack(geometry, c)
            end
        
            conv2 = abs(Icr - Icr_assumed) / Icr
            Icr_assumed = Icr
            c_assumed = "something as a function of Ωc"
        end

        Ie = getIe(Mcr, Mdec, M, Icr, Itr)
        δ_mid, δ_dev , _  = getDelta(Mat, Sec, f, Ie, Mi, em,fps)
        dps = dps0 - (δ_mid - δ_dev)

        fc = fps/Eps*c/(dps-c) + M/Itr*c  #concrete stress at top fiber
        fps_calc = getFps2(Mat, Sec, F , Ωc, c, dps, fc)
        conv1 = abs(fps_calc - fps_assumed) / fps_assumed
        fps = fps_calc
    end

    δ_mid, _ , e  = getDelta(Mat, Sec, f, Ie, Mi, em,fps)
        
    return δ_mid , fps, e
end



"""
currently not implemented, due to using ACI limit on fps instead.
"""
function ultimateRegime(M::Float64, fps_assumed::Float64, Icr_assumed::Float64, element::Element)
    @unpack Mat, Sec, F, type = element
    @unpack fc′, Ec, Eps, fpy = Mat
    @unpack Aps, Atr, Itr, Zb, em, es, Ls, Ld, L = Sec
    @unpack w, Mg, fr, r, fpe = F

    if Ld < Ls 
        K1 = Ls/L - 1
        K2 = Ls/Ls(Ld/L)^2 - (Ls/L)^2 
    elseif Ld ≥ Ls
        K1 = Ld/L - 1
        K2 = 0
    end

    if Sd/dps0 ≤ 15
        ks = 0.0096*(Sd/dps0)
    elseif Sd/dps0 > 15
        ks = 0.144
    end

    Ωu = 0.895- 1.364*(Ls/L) *dps0/h - ks

    fps = clamp(fpe + Ωu*Eps*ϵcu*(dps0/c-1) , 0, fpy)
    
    Φu = M/(Ec*Ie) 


    eu = (em + Φ *L^2/6 * (3 * Ld / L * (-K1) - 3 / 4 - K2)) / (1 - fps * Aps / (Ec * Itr) * (L^2 / 8 - L * Ld / 2 + Ld^2 / 2))

    dpsu = dps0 + Φ *L^2/6 * (3 * Ld / L * (-K1) - 3 / 4 - K2)  + fps*Aps/(Ec*Ie)*(eu*(L^2/8-L*Ld/2+Ld^2/2))
    #calculate until dps converge. 

    δmid = Φu*L^2/6*(3/4 -(L/Ls)^2) - fps*Aps/(Ec*Ie)*(eu*L^2/8 - (eu-es)*Ld^2/6)
    return δmid 
end