"""
Figure 6 Linear elastic uncracked regime.
***for the first guess of fps, fpe is preferred.
This function (workflow) includes.
(2) Tendon stress from the linear elastic uncracked section analysis.

if the tendon profile and the loaded points arent follow fig 4,\n 
equation (7) through (14) needs to be applied.

"""
function linearElasticUncrackRegime(M::Float64, fps_assumed::Float64, element::Element)::Float64
    @unpack Mat, Sec, F, type = element
    @unpack fc′, Ec, Eps, fpy = Mat
    @unpack Aps, Atr, Itr, Zb, em, es, Ls, Ld, L = Sec
    @unpack w, Mg, fr, r, fpe = F

    #(6)
    Ω = getΩ(Sec)

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
            top_e = em +M*L^2/(6*Ec*Itr)*(3*Ld/L*(-K1) - 3/4 - K2)
            bot_e = 1 - fps_assumed*Aps/(Ec*Itr)*(L^2/8 - L*Ld/2 + Ld^2/2)
            e =  top_e/bot_e 
        elseif type == "secondary" || type == "columns"
            k1 = 1/72 ; k2 = 1/72 ; 
            Δ = k1*M*L^2/(Ec*Itr) - k2*fps_assumed*Aps*e*L^2/(Ec*Itr)
            e = em-Δ
        else 
            println("Invalid section type")
            return 0.0
        end

        new_fps = fpe + Ω*M*e/(Itr*Ec/Eps + Aps*(r^2+e^2)*Ω)
        tol = abs(new_fps-fps_assumed)/fps_assumed
        fps_assumed = new_fps
    end

    fps = fps_assumed

    if type == "primary" 
        δ = M * L^2 / (6 * Ec * Itr) * (3 / 4 - (Ls / L)^2) - fps * Aps / (Ec * Itr) * (e * L^2 / 8 - (e - es) * Ls^2 / 6)
    elseif type == "secondary" || type == "columns"
        k3 = 1/8; 
        δ = 23*M*L^2/(216*Ec*Itr) - k3*fps*Aps*e*L^2/(Ec*Itr)
    end

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
function linearElasticCrackedRegime(M::Float64, fps_assumed::Float64, Icr_assumed::Float64, element::Element)::Float64
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
function ultimateRegime(M::Float64, fps_assumed::Float64, Icr_assumed::Float64, element::Element)::Float64
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





"""
(3) function getMcr(Mat::Material, Sec::Section, f::Loads, Ω::Float64)::Float64
Only calculate once
Find the moment limit for the liniear elastic uncracked section region.
"""
function getMcr(Mat::Material, Sec::Section, f::Loads, Ω::Float64)
    @unpack fc′, Ec, Eps, fpy = Mat
    @unpack Aps, Atr, Itr, Zb, em = Sec
    @unpack w, Mg, fr, r, fpe = f


    @show Mcre = Aps * fpe * (em + Zb / Atr) + (fr * Zb) # Cracking moment due to initial effective prestress (mcre)
    @assert Mcre > 0 "Error Mcre is less than 0."
    @show ΔMcr = (Aps * em * (em + Zb / Atr) * (Mcre - Mg)) / ((1 / Ω * Itr * Ec / Eps) + Aps * (r^2 - em * Zb / Atr)) # Moment due to stress increase in external tendons.
    @assert ΔMcr > 0 "Error ΔMce is less than 0."
    
    Mcr = Mcre + ΔMcr
    return Mcr
end

"""
(12)
"""
function getDelta(Mat::Material, Sec::Section, f::Loads, Itr::Float64, M::Float64, e::Float64, fps::Float64)

    @unpack fc′, Ec, fpy = Mat
    @unpack Aps, Atr, Zb, em, es, Ls, Ld, L = Sec
    @unpack w, mg, fr, r, fpe = f

    #displacement 
    # due to the PS force
    #at mid span
    δ_mid⁻ = fps * Aps / (Ec * Itr) * (em * L^2 / 8 - (em - es) * Ls^2 / 6)
    # at deviator 
    δ_dev⁻ = fps * Aps / (Ec * Itr) * (es * Ls^2 / 6 + em * (L * Ls / 2 - 2 / 3 * Ls^2))

    # due to the applied force
    #at the mid span
    δ_mid⁺ = M * L^2 / (6 * Ec * Itr) * (3 / 4 - (Ls / L)^2)
    # at a deviator
    δ_dev⁺ = M * L^2 / (6 * Ec * Itr) * (3 * (Ls / L) * (1 - Ls / L) - (Ls / L)^2)

    δ_mid = δ_mid⁺ - δ_mid⁻

    δ_mid_cal = M * L^2 / (6 * Ec * Itr) * (3 / 4 - (Ls / L)^2) - fps * Aps / (Ec * Itr) * (e * L^2 / 8 - (e - es) * Ls^2 / 6)
    # @show Itr
    # @show δ_mid, δ_mid_cal
    @assert abs(δ_mid - δ_mid_cal) < 1e-6
    δ_dev = δ_dev⁺ - δ_dev⁻

    Δ = δ_mid - (δ_dev⁺ - δ_dev⁻)
    @assert Δ == δ_mid - δ_dev

    K1 = Ls / L - 1
    K2 = 0.0
    Δcalc = M * L^2 / (6 * Ec * Itr) * (3 * (Ls / L) * K1 + 3 / 4 + K2) - fps * Aps * e / (Ec * Itr) * (L^2 / 8 - L * Ls / 2 + Ls^2 / 2)

    # @show Δ - Δcalc
    @assert abs(Δ - Δcalc) < 1e-9
    # e = (em + M*L^2/(6*Ec*Itr)*(3/4-(Ls/L)^2)) / (1 - fps*Aps/(Ec*Itr) * (L^2/8 - L*Ld/2 +Ld^2/2))
    e = (em + M * L^2 / (6 * Ec * Itr) * (3 * Ld / L * (-K1) - 3 / 4 - K2)) / (1 - fps * Aps / (Ec * Itr) * (L^2 / 8 - L * Ld / 2 + Ld^2 / 2))
    # println(e, e1)
    # @assert e == e1
    return δ_mid, δ_dev, e
end

"""
(19)
"""
function getFps1(Mat::Material, Sec::Section, f::Loads, Ωc::Float64, c::Float64, dps::Float64)

    @unpack fc′, Ec, Eps, fpy = Mat
    @unpack em, es, em0, dps0, Aps, Atr, Itr, Zb = Sec
    @unpack w, mg, fr, r, fpe, ϵpe, ϵce = f

    first_term = Eps * (ϵpe + Ωc * ϵce)
    second_term = Ωc * fc′ * Eps / Ec * (dps / c - 1)
    fps = first_term + second_term
    if fps > fpy
        return fpy
    else
        return fps
    end
end

"""
(23)
"""
function getFps2(Mat::Material, Sec::Section, f::Loads, Ωc::Float64, c::Float64, dps::Float64, fc::Float64)

    @unpack fc′, Ec, Eps, fpy = Mat
    @unpack em, es, em0, Aps, Atr, Itr, Zb = Sec
    @unpack w, mg, fr, r, fpe, ϵpe, ϵce = f
    #fc is suppose to be the stress in top concrete fiber, 
    # use constitution equation to get the stress in the top concrete fiber
    first_term = Eps * (ϵpe + Ωc * ϵce)
    second_term = Ωc * fc * Eps / Ec * (dps / c - 1)
    fps = first_term + second_term
    if fps > fpy
        return fpy
    else
        return fps
    end
end

"""
(20)
"""
function getLc(Sec::Section, Mcr::Float64, M::Float64)
    @unpack Ls, L = Sec
    Lc = L - 2 * Ls * Mcr / M
    if Lc < L - 2 * Ls
        return L - 2 * Ls
    else
        return Lc
    end
end

"""function getΩ(Sec::Section)
equation (6)
"""
function getΩ(Sec::Section)
    @unpack em, es, Ls, Ld, L = Sec

    if Ld < Ls 
        Ω = 1 - Ls/L + Ld^2*(es-em)/(3*L*Ls*em)
    elseif Ld ≥ Ls
        Ω = 1 - es/em *Ls/L + (es-em)/em*(Ls^3/(3*L*Ld)+Ld/L)
    end
    return Ω
end

"""
(21)
"""
function getΩc(Ω::Float64, Icr::Float64, Lc::Float64, Sec::Section)
    @unpack em, es, L, Ld, Ls, Itr = Sec

    if Ld < Ls
        if (L - 2 * Ls) < Lc < (L - 2 * Ld)
            Ωc = Ω * Icr / Itr + (1 - Icr / Itr) *
                                 (1 - L / (4 * Ls) + Lc / (2 * Ls) - Lc^2 / (4 * L * Ls) - Ls / L)
        elseif Lc >= L - 2 * Ld
            Ωc = Ω * Icr / Itr + (1 - Icr / Itr) *
                                 (1 - Ls / L - Ld^2 / (L * Ls) +
                                  (1 - es / em) * (L * Lc / (4 * Ld * Ls) - Lc^2 / (4 * Ld * Ls) +
                                                   Lc^3 / (12 * L * Ld * Ls) - L^2 / (12 * Ld * Ls) + 2 * Ld^2 / (3 * L * Ls)) +
                                  es / em * (Lc / (2 * Ls) - L / (4 * Ls) - Lc^2 / (4 * L * Ls) + Ld^2 / (L * Ls)))
        else
            println("Warning: Lc is out of range")
        end
    elseif Ld >= Ls
        Ωc = Ω * Icr / Itr + (1 - Icr / Itr) *
                             (1 - 2 * Ls / L +
                              (1 - es / em) * (L * Lc / (4 * Ld * Ls) - Lc^2 / (4 * Ld * Ls) +
                                               Lc^3 / (12 * L * Ld * Ls) - L^2 / (12 * Ld * Ls) + Ld / L - Ls^2 / (3 * L * Ld)) +
                              es / em * (Lc^2 / (4 * L * Ls) - L / (4 * Ls) + 2 * Ld / L - Ls / Ls))
    else
        println("Warning: Ld is out of range")
    end

    @assert Ωc > 0

    return Ωc
end

"""
(25)
"""
function getDeltamid()
    first_term = M * L^2 / (6 * Ec * Ie) * (3 / 4 - (Ls / L)^2)
    second_term = fps * Aps / (Ec * Ie) * (e * L^2 / 8 - (e - es) * Ld^2 / 6)
    return first_term + second_term
end

"""
(26)
Mdec is decompression moment
is the moment from externally post tension 
Mdec  = fpe*Aps*em
"""
function getIe(Mcr::Float64, Mdec::Float64, M::Float64, Icr::Float64, Itr::Float64)
    k = (Mcr - Mdec) / (M - Mdec)
    first_term = k^3 * Itr
    second_term = (1 - k^3) * Icr
    # @show first_term + second_term
    return clamp(abs(first_term + second_term), 0, Itr)
end

# """
# (27)
# This function is similar to linear uncrack
# """

"""
(28)
dps0 : initial effective post tension dendon depth
"""
function getDps(dps0::Float64, Δ::Float64)
    K1 = Ls / L - 1
    K2 = 0.0

    dps = dps0 + M * L^2 / (6 * Ec * Ie) * (3 * Ld / L * (-K1) - 3 / 4 - K2) +
          fps * Aps / (Ec * Ie) * e * (L^2 / 8 - L * Ld / 2 + Ld^2 / 2)
    @assert dps == dps0 + Δ
    return dps
end




"""
Fig 7 in the paper.
"""
function main()

    # I might have to unpack here.

    for i in eachindex(M)
        Mi = M[i]
        println(Mi)
        Lc = getLc(Sec, Mcr, Mi)
        loop1()

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

end

function loop1()

    conv1 = 1 #first convergence criteria
    counter1 = 0 #loop1 counter
    while conv1 > 1e-6
        counter1 += 1
        if counter1 > 1000
            println("Warning: 1st iteration did not converge")
            break
        end

        #assume value of Itr and fps
        # loop2()
        
        # println("Icr = ", Icr)
        # println("Ac_req ", Ac_req)
        # println("c: ", c)
        # @show Mcr , Mdec, Mi , Icr, Itr
        println(Mi)
        println(typeof(Mi))
        Ie = getIe(Mcr, Mdec, Mi, Icr, Itr)
        # println("Ie/Icr" , Ie/Icr)
        δ_mid, δ_dev, e = getDelta(Mat, Sec, f, Ie, Mi, em, fps)
        dps = dps0 - (δ_mid - δ_dev)
        fc = fps / Eps * c / (dps - c) + Mi / Itr * c
        # println("fc: ", fc)
        # @assert fc <= 0.003
        fps_calc = getFps2(Mat, Sec, f, Ωc, c, dps, fc)
        conv1 = abs(fps_calc - fps) / fps
        fps = fps_calc
        #plot convergence of fps, icr and dps using Makie
    end
end

    

function loop2()

    conv2 = 1 #second convergence criteria
    counter2 = 0 #loop2 counter
    while conv2 > 1e-6
        println("Im here in loop2")
        println(dps,fps,fpe,Ie,Itr)
        println(Ωc)
        println(Icr)
        # println("counter")
        counter2 += 1
        if counter2 > 1000
            println("Warning: 2nd iteration did not converge")
            break
        end
        Ωc = getΩc(Ω, Icr, Lc, Sec)
        # ps_force_i = Aps*fps

        #finding compression depth (c)
        c = getC()
        #calculate Icr
        Icr_calc = get_Icrack(c)

        conv2 = abs(Icr_calc - Icr) / Icr_calc

        Icr = Icr_calc
    end
end


function getC()
    c = 10.0 #dummy
    conv_c = 1
    counter_c = 0
    while conv_c > 1e-6
        counter_c += 1
        if counter_c > 1000
            println("Warning: 3rd iteration did not converge")
            break
        end
        #centroid of concrete area might not be at c/2
        Ac_req = Mi / (dps - c / 2) / (0.85 * fc′)

        # Ac_req = ps_force_i /0.85/fc′
        new_c = get_C(Ac_req)
        conv_c = abs(new_c - c) / new_c
        c = new_c
    end
end