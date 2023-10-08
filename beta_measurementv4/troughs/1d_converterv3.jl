include("/Users/alexjacoby/Documents/Research_Code/polarimeter/basic_functionality/1d_averager.jl")
import DelimitedFiles, Plots, LsqFit





index1 = DelimitedFiles.readdlm("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv4/troughs/Estimated_Angles.txt")
index2 = DelimitedFiles.readdlm("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv4/troughs/exact_angles.txt")

Plots.plot(index1[3, :]);
Plots.plot!(-index2[3, :])


small_fn_vec = vcat(reverse([string("-", i) for i in 1:20]), "0", [string("+", i) for i in 1:20])





cd("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv4/troughs/troughs_png/")

# n=11
# small_fn = small_fn_vec[n]

vl = length(small_fn_vec)
parms = zeros(vl, 6)
error = zeros(vl, 6)
data = zeros(2001, vl)

parms0 = [0.2, 0.1, 0.3, 2.5, 0, 0]

for (n, small_fn) in pairs(small_fn_vec)

    begin
        raw_image = Images.Gray.(Images.load(string(small_fn, ".png")))[270:710, 710:1340]
        dat = Float64.(raw_image)
        # Images.Gray.(correction*0.1)
        correction = (Float64.(Images.Gray.(Images.load("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv4/background_low_intensity.png"))[270:710, 710:1340])) .^ (-1)
        correction = (*(size(correction)...) / sum(correction)) * correction
        dat = dat .* correction
        # Images.Gray.(dat)
    end


    aangle = index2[2, :][n]
    (x, y) = convert_1d(dat, ϕ=-aangle)
    # x = x[250:1750]
    # y = y[250:1750]
    x = π * x / max(x...)


    Plots.plot(x, y)

    model(r::Vector{Float64}, parms::Vector{Float64}) = parms[1] * cos.(parms[4] * r .- parms[5]) .+ parms[2] * cos.(2 * parms[4] * r .+ parms[6] .- 2 * parms[5]) .+ parms[3]
    fit = LsqFit.curve_fit(model, x, y, parms0)

    parms[n, :] = fit.param
    error[n,:] = LsqFit.estimate_errors(fit)
    parms0 = fit.param

    data[:, n] = y
    # Plots.plot(x, y);
    # Plots.plot!(x, model(x, parms0));
    # Plots.plot!(x, model(x, fit.param))
end





model(r::Vector{Float64}, parms::Vector{Float64}) = parms[1] * cos.(parms[4] * r .- parms[5]) .+ parms[2] * cos.(2 * parms[4] * r .+ parms[6] .- 2 * parms[5]) .+ parms[3]
x = Vector(-π:(2π/(2000)):π)
ani = Plots.@animate for n in 1:vl
    Plots.plot(x, data[:, n], seriestype=:scatter, legend=:topright, ms=2.0, color=:black, markerstrokecolor=:match, label = "Data", title = "Trough Matched Region")
    Plots.plot!(x, model(x, parms[n, :]), linewidth=5, label = "Fitted")
end
Plots.gif(ani, "/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv4/troughs/Trough_Fullrange.gif")




Plots.plot(index2[3, :], parms[:, 6] .+ 3* 2π, yerror = error[:,6], seriestype = :scatter, label = "Relative Phase", title = "Relative Phases for Trough Matched Point")
linear_model(x::Vector{Float64}, parms2::Vector{Float64}) = parms2[2] .+ parms2[1] * x
p0 = [0.0, 0]
linear_fit = LsqFit.curve_fit(linear_model, index2[3, :], parms[:, 6] .+ 3* 2π, p0)
p = linear_fit.param
-p[2] / p[1]
Plots.plot!(index2[3,:], linear_model(index2[3,:], p), linewidth = 4, label = "Fitted")


Plots.savefig("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv4/troughs/Zero_Crossing.pdf")