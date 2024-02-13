#visualize the final catalog. 
using Makie, GLMakie

results = CSV.read("src/Catalogs/FEB6_4_catalog_static.csv", DataFrame); 

f_catalog = Figure(size = (30000,20000))
titles =  ["fc′", "as", "dps", "fpe", "carbon", "L", "t", "Lc", "T", "Pu", "Mu", "Vu"]
Axes = Vector{Axis}(undef, length(titles))
color_map = x-> x==0 ? :black : :green


#special filter
# results = subset(results , :fc′ => x-> x .== 60)

results_no_zeros = subset(results, :Mu => x-> x .!= 0)
results_zeros = subset(results, :Mu => x-> x  .== 0)
results_zeros = subset(results_zeros, :dps => x-> x .!= 0)

for i in eachindex(titles)
    row = div(i-1,5)+1
    col= mod(i+4,5)+1
    # println(row,",", col)
    Axes[i] = Axis(f_catalog[row, col], title = titles[i], xlabel = "dps [mm]")
    #plot non-zero as colorful plot
    scatter!(Axes[i], results_no_zeros[!, :dps], results_no_zeros[!, titles[i]], color = results_no_zeros[!,:Mu] , colormap = :plasma, markersize = 7.5)
    #plot zeros as black on top (that's why it is separated and plotted later.)
    scatter!(Axes[i], results_zeros[!, :dps], results_zeros[!, titles[i]], color = :black , marker = 'x', markersize = 20)
end
