using Makie, GLMakie
using CSV, DataFrames
using Formatting
using Statistics

save_file = true;

#load the dataset. (could be anything.)
path = "src//Tables//Compiled Concrete EPD Data Revised 11_06_2023.csv"
dataset_Broyles = CSV.read(path, DataFrame);

# Extract fc′ and embodied carbon datapoints from the loaded DataFrame
fc′_column_name = "Concrete Compressive Strength (psi)"
ec_column_name = "A1-A3 Global Warming Potential (kg CO2-eq)"

fc′_Broyles_psi = dataset_Broyles[!, fc′_column_name    ]
ec_Broyles = dataset_Broyles[!, ec_column_name]

#number formatting for clean visualization
printx(x,N) = sprintf1("%10.$(N)f",x)

f1 = Figure(size = (800,800))
ax1 = Axis(f1[1,1],
title = "Broyles'",
xlabel = "Concrete Strength [MPa]",
ylabel = "Embodied Carbon [kg CO2e/m3]",
xminorticks = 0:2500:20000,
xminorgridvisible = true,
yminorticks = 0:500:4000,
yminorgridvisible = true,
xtickformat =  x -> printx.(x,0),
titlesize = 30,
xlabelsize = 20,
ylabelsize = 20,
)

scatter!(ax1, fc′_Broyles_psi, ec_Broyles)
text!(ax1, string(length(fc′_Broyles_psi)) *" Datapoints", position=(10000,3000), fontsize = 25)
f1

if save_file
save("2002_Overall_plot.png", f1)
end

# Note the psi and MPa conversion
#convert to MPa [1 psi = 0.00689476 MPa]
fc′_Broyles_MPa = fc′_Broyles_psi*0.00689476;



#define groups as done in the paper
label = Vector{Int64}(undef, length(fc′_Broyles_MPa))
for i in eachindex(fc′_Broyles_MPa)
    fc′ = fc′_Broyles_psi[i]
    if fc′ <=2000
        label[i] = 1
    elseif fc′<= 3000
        label[i] = 2
    elseif fc′<= 4000
        label[i] = 3
    elseif fc′<=5000
        label[i] = 4
    elseif fc′<=6000
        label[i] = 5
    elseif fc′<=8000
        label[i] = 6
    elseif fc′<= 10000
        label[i] = 7
    elseif fc′<= 12000
        label[i] = 8
    else
        label[i] = 9
    end
end

# add the label into the dataset. 
dataset_Broyles[!,"Label"] = label;

n_label = combine((groupby(dataset_Broyles, :Label, sort=true)), nrow)[!, :nrow]

#define text position for labeling
text_pos = Point2f.([[i ,100.] for i in 1.:9.])


f3 = Figure(size = (1000,400))
ax3 = Axis(f3[1,1],
title = "fc′ vs Embodied Carbon",
titlesize = 20,
xticks=  1:9,
xlabel = "Concrete Strength [psi]",
xminorgridvisible = true,
xtickformat =  x -> ["< 2,000",
"2,001-3,000",
"3,001-4,000",
"4,001-5,000",
"5,001-6,000",
"6,001-8,000",
"8,001-10,000",
"10,001-12,000",
">12,000"][Int.(x)],
ylabel = "Embodied Carbon [kgCO2e/m3]",
yticks = 0:100:800)
boxplot!(ax3,label,ec_Broyles, color = :red, show_outliers = false)
text!(ax3, string.(n_label), position= text_pos)
f3

#get company names
company_names = unique(dataset_Broyles[!, :Company])
company_n_EPDs = Dict{String, Int64}()
n_company = length(company_names)
for i in company_names
    company_n_EPDs[i] = size(dataset_Broyles[dataset_Broyles.Company .== i,:],1)
end
number_of_EPDs = [ company_n_EPDs[i] for i in company_names]
perm = sortperm(number_of_EPDs, rev= true)
company_names = company_names[perm]

f_company = Figure(size = (1000,1000))
a_company = Axis(f_company[1,1], xticks = (1:n_company, company_names), xticklabelrotation = 45, 
title = "Number of EPDs per company",
titlesize = 20, 
xlabel = "Company",
ylabel = "Number of EPDs")
barplot!(a_company,
   1:n_company,         
   number_of_EPDs[perm],
   bar_labels = number_of_EPDs[perm],
   gap = 0.5)
f_company



#Plot trend fc′ and EC

max_y = maximum(ec_Broyles)
text_pos = Point2f.([[i ,100.] for i in 1.:9.])
fc′_design = 0.00689476*[2000.0, 2500.0, 3500.0, 4500.0, 5500.0, 7000.0, 9000.0, 11000, 12000.0] #MPa

f4 = Figure(size = (1200,500))
ax4 = Axis(f4[1,1],
title = "EC trend schemes", titlesize = 20,
xlabel = "Concrete Strength [MPa]",
ylabel = "Embodied Carbon [kgCO2e/m3]",

xlabelsize = 15,
ylabelsize = 15,

xticklabelrotation = 45,

limits = (0, nothing, 0, 1500),

xticks = fc′_design, xtickformat =  x -> printx.(x,2),
xminorgridvisible = true,
yticks = 0:100:max_y)

right_grid = f4[1,2] = GridLayout()
menu = Menu(right_grid[1,1], options = company_names, tellheight = false , tellwidth = false)

line1 = Vector{Float64}(undef, 9)
line2 = Vector{Float64}(undef, 9)
line3 = Vector{Float64}(undef, 9)
#plot static part of the figure.
for i in 1:9
    sub_dataset = dataset_Broyles[dataset_Broyles.Label .== i, ec_column_name ]
    first_quantile  = quantile!(sub_dataset, 0.25)
    second_quantile = quantile!(sub_dataset, 0.50)
    third_quantile  = quantile!(sub_dataset, 0.75)
    IQR = third_quantile - first_quantile
    line1[i] = second_quantile
    line2[i] = third_quantile + 1.5*IQR
    line3[i] = clamp(first_quantile - 1.5*IQR, minimum(sub_dataset), Inf)
end 

scatter!(ax4, fc′_Broyles_MPa, ec_Broyles, color =color = colorant"#D3D3D3", label = "All datapoints")
lines!(ax4, fc′_design,line2, label = "Maximum", linewidth = 3, color = :red)
lines!(ax4, fc′_design,line1, label = "Median",  linewidth = 3, color = :green)
lines!(ax4, fc′_design,line3, label = "Minimum", linewidth = 3, color = :blue)
Legend(right_grid[2,1],ax4)

line4 = Observable(Vector{Float64}(undef, 9));
points_x = Observable(Vector{Float64}())
points_y = Observable(Vector{Float64}())
#select a company
company_name = on(menu.selection) do name
    company_dataset = dataset_Broyles[dataset_Broyles.Company .== name, :]
    points_x.val = company_dataset[!, fc′_column_name]*0.00689476
    points_y[] = company_dataset[!, ec_column_name]
    
    dummy_line4 = Vector{Float64}(undef, 9);

    for i in 1:9
        company_dataset_label = company_dataset[company_dataset.Label .== i, ec_column_name ]
        if length(company_dataset_label) > 0
            @show val = median(company_dataset[company_dataset.Label .== i, ec_column_name ])
            if val == 0 
                dummy_line4[i] = nothing 
            else 
                dummy_line4[i] = val
            end
        else
            println(i)
        end
    end
    line4[] = dummy_line4 
end


#plot the previous data, greyed out
scatter!(ax4, points_x, points_y, color = :black)
scatter!(ax4, fc′_design,line4, label = company_name, marker = :utriangle , markersize = 20, color = :orange)
lines!(ax4, fc′_design, line4, color = :skyblue)
# Legend(right_grid[2,1],ax4)
f4

#now, apply beam calculation on median points AND scatter points.

# include("reinforced_concrete.jl")
# include("beam_optimization.jl")



