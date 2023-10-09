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






import Images, FileIO, LinearAlgebra, LsqFit, ImageTransformations, OffsetArrays, Meshes, Plots
include("/Users/alexjacoby/Documents/Research_Code/polarimeter/basic_functionality/Rotation_Rectangle.jl")
include("/Users/alexjacoby/Documents/Research_Code/polarimeter/basic_functionality/1d_averager.jl")





n = 6
small_fn = small_fn_vec[n]
index = DelimitedFiles.readdlm("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv4/troughs/exact_angles.txt")
γ = -31.5 * (π / 180)
θ = index[2, n] - γ
θ_deg = θ * (180 / π)
begin
    raw_imagepr = Images.Gray.(Images.load(string("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv4/troughs/troughs_png/", small_fn, ".png")))

    raw_image = OffsetArrays.centered(ImageTransformations.imrotate(raw_imagepr, γ))[-450:50, -70:530]
    dat = Float64.(raw_image)

    correction_imagepr = Images.Gray.(Images.load("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv4/background_low_intensity.png"))
    correction_image = OffsetArrays.centered(ImageTransformations.imrotate(correction_imagepr, γ))[-450:50, -70:530];

    correction = Float64.(correction_image) .^ (-1)
    correction = (*(size(correction)...) / sum(correction)) * correction;
    dat = dat .* correction
    image = Images.Gray.(dat);
end






(N, M) = size(image)
LX = M - 1
LY = N - 1

Aspect = LX / LY
(R, S) = rotation_crop(dims=size(image), angle=θ, AR=Aspect)
rotated = OffsetArrays.centered(ImageTransformations.imrotate(image, θ));
# Plots.plot(rotated);
# Plots.plot!([S, -S, S, -S], [R, R, -R, -R], seriestype=:scatter, color=:red)


sample = Float64.(rotated[-R:R, -S:S])
Images.Gray.(sample)

strip_length = 30
strip_skip = 1

partitions = ((size(sample)[1] - strip_length) ÷ strip_skip) + 1
partition_rngs = [(strip_skip*(n-1)+1):(strip_skip*(n-1)+strip_length) for n in 1:partitions]

partition_rngs[length(partition_rngs)]
data_vec = [dropdims(sum(sample[rng, :], dims=(1)), dims=(1)) / length(rng) for rng in partition_rngs]





model(r::Vector{Float64}, parms::Vector{Float64}) = parms[1] * cos.(parms[4] * r .- parms[5]) .+ parms[2] * cos.(2 * parms[4] * r .+ parms[6] .- 2 * parms[5]) .+ parms[3]



params = zeros(partitions, 6)


parms0 = [0.1, 0.05, 0.2, 2.0, 0, 0]
for n in 1:partitions
    y = data_vec[n]
    x = π * Vector(-1.0:2/(length(y)-1):1)

    fit = LsqFit.curve_fit(model, x, y, parms0)
    params[n, :] = fit.param
    parms0 = fit.param
    # Plots.plot(x, y, label="Data", seriestype=:scatter, legend=:topright)
    # Plots.plot!(x, model(x, fit.param), label="Fitted", linewidth=5)
end
@. linear_model(x, parms2) = parms2[2] + parms2[1] * x
p0 = [0.0, 0]
linear_fit = LsqFit.curve_fit(linear_model, Vector{Float64}(1:partitions), params[:, 5], p0)
slope = linear_fit.param[1]
error = LsqFit.estimate_errors(linear_fit)[1]
Plots.plot(1:partitions, params[:, 5]);
Plots.plot!(1:partitions, linear_model(Vector(1:partitions), linear_fit.param))

Plots.plot(1:partitions, params[:, 1])

animation = Plots.@animate for n in 1:partitions
    y = data_vec[n]
    x = π * Vector(-1.0:2/(length(y)-1):1)
    parms0 = params[n, :]
    fit = LsqFit.curve_fit(model, x, y, parms0)
    Plots.plot(x, y, label="Data", seriestype=:scatter, legend=:topright)
    Plots.plot!(x, model(x, fit.param), label="Fitted", linewidth=5)
end


Plots.gif(animation, "moving.gif", fps=15)






























θ0 = 0 * (pi / 180)
range = 2 * (pi / 180)
increments = 40
slope_vec1 = []
error_vec1 = []
θ_vec = Vector((θ0-range):(2*range/(increments-1)):(θ0+range))
for θ in θ_vec
    (N, M) = size(image)
    LX = M - 1
    LY = N - 1
    Aspect = LX / LY
    (R, S) = rotation_crop(dims=size(image), angle=θ, AR=Aspect)
    rotated = OffsetArrays.centered(ImageTransformations.imrotate(image, θ))
    # Plots.plot(rotated);
    # Plots.plot!([S, -S, S, -S], [R, R, -R, -R], seriestype=:scatter, color=:red)


    sample = Float64.(rotated[-R:R, -S:S])



    strip_length = 30
    strip_skip = 1

    partitions = ((size(sample)[1] - strip_length) ÷ strip_skip) + 1
    partition_rngs = [(strip_skip*(n-1)+1):(strip_skip*(n-1)+strip_length) for n in 1:partitions]

    partition_rngs[length(partition_rngs)]
    data_vec = [dropdims(sum(sample[rng, :], dims=(1)), dims=(1)) / length(rng) for rng in partition_rngs]


    model(r::Vector{Float64}, parms::Vector{Float64}) = parms[1] * cos.(parms[4] * r .- parms[5]) .+ parms[2] * cos.(2 * parms[4] * r .+ parms[6] .- 2 * parms[5]) .+ parms[3]

    params = zeros(partitions, 6)
    parms0 = [0.1, 0.05, 0.2, 2.0, 0, 0]
    for n in 1:partitions
        y = data_vec[n]
        x = π * Vector(-1.0:2/(length(y)-1):1)

        fit = LsqFit.curve_fit(model, x, y, parms0)
        params[n, :] = fit.param
        parms0 = fit.param
        # Plots.plot(x, y, label="Data", seriestype=:scatter, legend=:topright)
        # Plots.plot!(x, model(x, fit.param), label="Fitted", linewidth=5)
    end

    linear_model(x::Vector{Float64}, parms2::Vector{Float64}) = parms2[2] .+ parms2[1] * x
    p0 = [0.0, 0]
    linear_fit1 = LsqFit.curve_fit(linear_model, Vector{Float64}(1:partitions), params[:, 5], p0)
    slope1 = linear_fit1.param[1]
    error1 = LsqFit.estimate_errors(linear_fit1)[1]


    append!(slope_vec1, slope1)
    append!(error_vec1, error1)

    println(string("θ = ", θ, " done"))
end




Plots.plot((180 / pi) * (θ_vec .+ γ), slope_vec1, yerror=error_vec1, legend=:none, xlabel="Angle (degrees)", ylabel="Phase Slope (Arbitrary Scale)", seriestype=:scatter)

θ_vec
slope_vec1
linear_model(x::Vector{Float64}, parms2::Vector{Float64}) = parms2[2] .+ parms2[1] * x
p0 = [0.0, 0]
final_fit = LsqFit.curve_fit(linear_model, θ_vec, slope_vec1, p0)
final_parms = final_fit.param
finalangle = -final_parms[2] / final_parms[1]
(finalangle + γ) * (180 / pi)

Plots.plot((180 / pi) * θ_vec, slope_vec1, yerror=error_vec1, legend=:none, xlabel="Angle (degrees)", ylabel="Phase Slope (Arbitrary Scale)", seriestype=:scatter);
Plots.plot!((180 / pi) * θ_vec, linear_model(θ_vec, final_parms))


