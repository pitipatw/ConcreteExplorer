#visualize the final catalog. 
using Makie, GLMakie, CairoMakie
using CSV, DataFrames

function catalog_vis(catalog_path::String)
    GLMakie.activate!()
    println("Visualizing a catalog from $catalog_path...")
    #load the catalog.
    # results = CSV.read("src/Catalogs/20_03_catalog_for_paper.csv", DataFrame);
    results = CSV.read(catalog_path, DataFrame)
    color_range = (minimum(results[!, :Mu]), maximum(results[!, :Mu]))
    x_axis = "fc′"
    x_label = "fc′ [MPa]"
    x_min = minimum(results[!, x_axis])
    x_max = maximum(results[!, x_axis])

    f_catalog = Figure(size=(3000, 2000))

    titles = ["fc′", "as", "dps", "fpe", "fps", "dosage", "carbon", "L", "t", "Lc", "T", "Pu", "Mu", "Vu"]
    units = ["MPa", "mm2", "mm", "MPa", "MPa", "kg steel/m3 concrete", "kgCO2e/m", "mm", "mm", "mm", "-", "kN", "kNm", "kN"]

    Axes = Vector{Axis}(undef, length(titles))
    g1 = f_catalog[3, 4] = GridLayout()
    g2 = f_catalog[3, 5] = GridLayout()

    slides = Vector{Any}(undef, length(titles))
    for i in eachindex(titles)
        unit = units[i]
        min = minimum(results[!, titles[i]])
        max = maximum(results[!, titles[i]])
        gap = (max - min) / 10
        slide_row = div(i - 1, 5) + 1
        slide_col = 2 * (mod(i + 4, 5) + 1)
        if gap == 0
            range = min
        else 
            range = LinRange(min, max, 1000)
        end
        #Box is for checking the plot area.
        # Box(f_catalog[slide_row, slide_col], color = (:red, 0.2), strokewidth = 0)

        slides[i] = IntervalSlider(f_catalog[slide_row, slide_col], range=range, startvalues=(min, max),
            horizontal=false) #, tellheight =false)

        # if i <= 6
        #     slides[i] = IntervalSlider(g1[i,1], range=range, startvalues= (min,max), tellheight =false)
        # elseif i > 6
        #     slides[i] = IntervalSlider(g2[i-6,1], range=range, startvalues= (min,max), tellheight =false)
        # end
    end

    #based, plot everything in grey first.
    for i in eachindex(titles)
        row = div(i - 1, 5) + 1
        col = mod(i + 4, 5) + 1

        #Plot the graphs on odd rows
        col = 2 * (col - 1) + 1
        # println(row,",", col)
        # Box is here for checking the plot area.
        # Box(f_catalog[row, col], color = (:red, 0.2), strokewidth = 0)
        Axes[i] = Axis(f_catalog[row, col], title=titles[i], xlabel=x_label, xticks=0:25:350, limits=(-10, x_max * 1.05, nothing, nothing))
        #plot non-zero as colorful plot
        scatter!(Axes[i], results[!, x_axis], results[!, titles[i]], color=:grey, markersize=2, transparent=true)
        #plot zeros as black on top (that's why it is separated and plotted later.)
    end

    #special filter
    # results = subset(results , :fc′ => x-> x .== 60)
    # results = subset(results , :fc′ => x-> x .<= 40)
    # results = subset(results , :as => x-> x .<= 200)

    # results_no_zeros = subset(results, :Mu => x -> x .!= 0)
    # results_zeros = subset(results, :Mu => x -> x .== 0)
    # results_zeros = subset(results_zeros, :dps => x -> x .!= 0)


    val = Vector{Any}(undef, 14)
    for i in eachindex(val)
        val[i] = slides[i].interval
    end

    dftovec(df::DataFrame) = [collect(df[i, titles]) for i in 1:size(df)[1]]

    points = dftovec(results)

    colors = lift(val...) do fc′, as, dps, fpe, fps, dosage, carbon, L, t, Lc, T, Pu, Mu, Vu
        map(points) do p
            x = (fc′[1] <= p[1] <= fc′[2]) &&
                (as[1] <= p[2] <= as[2]) &&
                (dps[1] <= p[3] <= dps[2]) &&
                (fpe[1] <= p[4] <= fpe[2]) &&
                (fps[1] <= p[5] <= fps[2]) &&
                (dosage[1] <= p[6] <= dosage[2]) &&
                (carbon[1] <= p[7] <= carbon[2]) &&
                (L[1] <= p[8] <= L[2]) &&
                (t[1] <= p[9] <= t[2]) &&
                (Lc[1] <= p[10] <= Lc[2]) &&
                (T[1] <= p[11] <= T[2]) &&
                (Pu[1] <= p[12] <= Pu[2]) &&
                (Mu[1] <= p[13] <= Mu[2]) &&
                (Vu[1] <= p[14] <= Vu[2])
        end
    end

    marker_size = lift(colors) do c
        c .* 7.5
    end

    # colors = lift(colors) do c 
    #     c.*results[!, :Mu]
    # end

    for i in eachindex(titles)
        row = div(i - 1, 5) + 1
        col = mod(i + 4, 5) + 1
        col = 2 * (col - 1) + 1
        # println(row,",", col)
        # Axes[i] = Axis(f_catalog[row, col], title=titles[i], xlabel="dps [mm]", xticks=0:25:350, limits=(0, 310, nothing, nothing))
        #plot non-zero as colorful plot
        scatter!(Axes[i], results[!, x_axis], results[!, titles[i]], colormap=:plasma, markersize=marker_size, colorrange=color_range)
        #plot zeros as black on top (that's why it is separated and plotted later.)
        # scatter!(Axes[i], results_zeros[!, :dps], results_zeros[!, titles[i]], color=:black, marker='x', markersize=20)
    end

    # for i in eachindex(titles)
    #     row = div(i - 1, 5) + 1
    #     col = mod(i + 4, 5) + 1
    #     # println(row,",", col)
    #     # Axes[i] = Axis(f_catalog[row, col], title=titles[i], xlabel="dps [mm]", xticks=0:25:350, limits=(0, 310, nothing, nothing))
    #     #plot non-zero as colorful plot
    #     scatter!(Axes[i], results_no_zeros[!, :dps], results_no_zeros[!, titles[i]], color=results_no_zeros[!, :Mu], colormap=:plasma, markersize=7.5)
    #     #plot zeros as black on top (that's why it is separated and plotted later.)
    #     scatter!(Axes[i], results_zeros[!, :dps], results_zeros[!, titles[i]], color=:black, marker='x', markersize=20)
    # end

    return f_catalog
end

# catalog_vis("src//Catalogs//18_03_catalog_for_test.csv")