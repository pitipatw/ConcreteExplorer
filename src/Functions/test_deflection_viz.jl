include("get_Deflection.jl")

element1 = PixelFrameElement()
outputs = get_Deflection(element1, 4.448*8300)

fps_history = outputs["fps_history"]
dps_history = outputs["dps_history"]
Icr_history = outputs["Icr_history"]
Ie_history  = outputs["Ie_history" ]
c_history   = outputs["c_history"]
dis_history = outputs["dis_history"]
fc_history  = outputs["fc_history"]
dis_dev_history= outputs["dis_dev_history"]

P = outputs["P"]
Mcr =outputs["Mcr"]
Mdec = outputs["Mdec"]
Ls = 502.7


#put this into another functin to be compared with.
df = CSV.File(joinpath(@__DIR__,"pixelframe_beam1.csv"))
df = DataFrame(df)
test_P = df[!,2]
test_d = df[!,3]

using Makie, GLMakie 
#Vizualization
fig1 = Figure(backgroundcolor = RGBf(0.98,0.98,0.98) , size = (1500, 1500))
fig2 = Figure(backgroundcolor = RGBf(0.98,0.98,0.98) , size = (2000, 1500))
ga = fig1[1,1] = GridLayout()
# gb = fig1[1,2] = GridLayout()
gb = fig2[1,1] = GridLayout()
title_name1 = [ "dps", "fps", "DisMid", "c", "Inertia(Crack(blue), eff(red))"] 
title_name2 = [ "Original", "Shifted by Moment crack"]
fig_monitor = Figure(resolution = (1200, 2000))
x1 = ["P [N]", "P [N]", "P [N]", "P [N]", "P [N]"]
y1 = ["dps [mm]", "fps [MPa]", "DisMid [mm]", "c [mm]", "Inertia [mm4]"]
x2 = ["Displacement [mm]", "Displacement [mm]"]
y2 = ["P [N]", "P [N]"]
axis_monitor1 = [Axis(ga[i,1], title = title_name1[i],ylabel = y1[i], xlabel = x1[i]) for i in 1:5]
axis_monitor2 = [Axis(gb[i,1], title = title_name2[i],ylabel = y2[i], xlabel = x2[i], yticks = -400000.:2500:40000)  for i in 1:2]

scatter!(axis_monitor1[1], P, dps_history, color = :red)
scatter!(axis_monitor1[2], P, fps_history, color = :red)
scatter!(axis_monitor1[2], P, fc_history, color = :blue)
scatter!(axis_monitor1[3], P, dis_history, color = :red)
scatter!(axis_monitor1[3], P, dis_dev_history, color = :blue)
scatter!(axis_monitor1[4], P, c_history, color = :red)
scatter!(axis_monitor1[5], P, Ie_history, color = :red ,label = "Ie")
scatter!(axis_monitor1[5], P, Icr_history, color = :blue, label= "Icr")
#add verticle line on each plot for Mcr
for i in 1:5
    vlines!(axis_monitor1[i], [Mcr*2/Ls], color = :black, label = "Mcr", linewidth = 5)
    #add verticle line for Mdec
    vlines!(axis_monitor1[i], [Mdec*2/Ls], color = :green, label = "Mdec")
end

# convert to in to mm
test_d = test_d .* 25.4
test_P = test_P .* 4.44822

plot!(axis_monitor2[1],dis_history[1:end],P[1:end], label = "calc", color = :blue)
plot!(axis_monitor2[1],test_d,test_P, label = "test", color = :red)

plot!(axis_monitor2[2],dis_history[1:end].-0.4,P[1:end].-Mcr*2/Ls, label = "calc", color = :blue)
plot!(axis_monitor2[2],test_d,test_P, label = "test", color = :red)
# display(fig_monitor)


fig3 = Figure(size = (800, 600))
ax3 = Axis(fig3[1, 1], ylabel = "Force Diff [N]", xlabel = "Displacement [mm]")
plot!(ax3,dis_history[1:end],P[1:end].-Mcr*2/Ls, label = "calc", color = :blue)
plot!(axis_monitor2[2], dis_history, dis_history.*1000, label = "dis/1000", color = :green, markersize = 1)
#plot 


display(fig1)
display(fig2)
# display(fig3)