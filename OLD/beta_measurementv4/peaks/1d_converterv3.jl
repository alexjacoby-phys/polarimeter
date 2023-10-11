import Images, FileIO, LinearAlgebra, LsqFit, ImageTransformations, OffsetArrays, Meshes, Plots, DelimitedFiles

include("/Users/alexjacoby/Documents/Research_Code/polarimeter/basic_functionality/Rotation_Rectangle.jl")
include("/Users/alexjacoby/Documents/Research_Code/polarimeter/basic_functionality/1d_averager.jl")


index = DelimitedFiles.readdlm("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv4/peaks/exact_angles.txt")


cd("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv4/peaks/peaks_png/")




small_fn_vec = vcat(reverse([string("-", i) for i in 1:30]), "0", [string("+", i) for i in 1:20])





n=23
small_fn = small_fn_vec[n]





γ = 12.6 * (π / 180)
begin
    raw_imagepr = Images.Gray.(Images.load("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv4/peaks/peaks_png/-9.png"));

    raw_image = OffsetArrays.centered(ImageTransformations.imrotate(raw_imagepr, γ))[-240:280, 20:600];
    dat = Float64.(raw_image)

    correction_imagepr = Images.Gray.(Images.load("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv4/background_low_intensity.png"))
    correction_image = OffsetArrays.centered(ImageTransformations.imrotate(correction_imagepr, γ))[-240:280, 20:600]

    correction = Float64.(correction_image) .^(-1)
    correction = (*(size(correction)...) / sum(correction)) * correction;
    dat = dat .* correction
    image = Images.Gray.(dat)
end


θ = (index[2,n] - γ)
(N, M) = size(image)
LX = M - 1
LY = N - 1
#Aspect = 1.2
#Aspect = 2.5
Aspect = LX / LY
(R, S) = rotation_crop(dims=size(image), angle=θ, AR=Aspect)
rotated = OffsetArrays.centered(ImageTransformations.imrotate(image, θ));
# Plots.plot(rotated);
# Plots.plot!([S, -S, S, -S], [R, R, -R, -R], seriestype=:scatter, color=:red)


sample = Float64.(rotated[-R:R, -S:S])
Images.Gray.(sample)


strip_length = 10
strip_skip = 2

partitions = ((size(sample)[1] - strip_length) ÷ strip_skip) + 1
partition_rngs = [(strip_skip*(n-1)+1):(strip_skip*(n-1)+strip_length) for n in 1:partitions]

partition_rngs[length(partition_rngs)]
data_vec = [dropdims(sum(sample[rng, :], dims=(1)), dims=(1)) / length(rng) for rng in partition_rngs]





model(r::Vector{Float64}, parms::Vector{Float64}) = parms[1] * cos.(parms[4] * r .- parms[5]) .+ parms[2] * cos.(2 * parms[4] * r .+ parms[6] .- 2 * parms[5]) .+ parms[3]



params = zeros(partitions, 6)


parms0 = [0.1, 0.05, 0.2, 2.5, 0, 0]
for n in 1:partitions
    y = data_vec[n]
    x = π * Vector(-1.0:2/(length(y)-1):1)
    parms0 = [0.1, 0.05, 0.2, 2.5, 0, 0]

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




animation = Plots.@animate for n in 1:partitions
    parms0 = params[n, :]
    y = data_vec[n]
    x = π * Vector(-1.0:2/(length(y)-1):1)
    fit = LsqFit.curve_fit(model, x, y, parms0)
    params[n, :] = fit.param
    Plots.plot(x, y, label="Data", seriestype=:scatter, legend=:topright)
    Plots.plot!(x, model(x, fit.param), label="Fitted", linewidth=5)
end

Plots.gif(animation, "moving.gif", fps=15)

params
params6_edited = [ cf_2f(params[n,:])[6] for n in 1:partitions]
Plots.plot(1:partitions, params6_edited)