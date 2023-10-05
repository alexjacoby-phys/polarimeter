import Images, FileIO, LinearAlgebra, LsqFit, ImageTransformations, OffsetArrays, Meshes
include("/Users/alexjacoby/Documents/Research_Code/polarimeter/basic_functionality/Rotation_Rectangle.jl")


















raw_image = Images.Gray.(Images.load("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv4/troughs/troughs_png/-5.png"))[270:710, 700:1350];

dat = Float64.(raw_image)

correction = (Float64.(Images.Gray.(Images.load("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv4/background_low_intensity.png"))[270:710, 700:1350])) .^ (-1)
correction = (*(size(correction)...) / sum(correction)) * correction;
dat = dat .* correction
Images.Gray.(dat)
function rotation(angle::Float64)
    return [cos(angle) -sin(angle); sin(angle) cos(angle)]
end










θ = -32.6*(π/180)
(N, M) = size(raw_image)
LX = M - 1
LY = N - 1
#Aspect = 1.2
#Aspect = 2.5
Aspect = 2.5
(R, S) = rotation_crop(dims=size(raw_image), angle=θ, AR=Aspect)


rotated = OffsetArrays.centered(ImageTransformations.imrotate(raw_image, θ))
Plots.plot(rotated);
Plots.plot!([S, -S, S, -S], [R, R, -R, -R], seriestype=:scatter, color=:red)


sample = Float64.(rotated[-R:R, -S:S])



strip_length = 10
strip_skip = 1

partitions = ((size(sample)[1] - strip_length) ÷ strip_skip) + 1
partition_rngs = [(strip_skip*(n-1)+1):(strip_skip*(n-1)+strip_length) for n in 1:partitions]

partition_rngs[length(partition_rngs)]
data_vec = [dropdims(sum(sample[rng, :], dims=(1)), dims=(1)) / length(rng) for rng in partition_rngs]





model(r::Vector{Float64}, parms::Vector{Float64}) = parms[1] * cos.(parms[4] * r .+ parms[5]) .+ parms[2] * cos.(2 * parms[4] * r .+ parms[6]) .+ parms[3]





n=1






y = data_vec[n]
x = π*Vector(-1.:2/(length(y)-1):1)
parms0 = [0.04, 0.06, 0.31, 2.37, 0, 0]# 0.12899908028433657-π, 0.48438720137510083 ]

fit = LsqFit.curve_fit(model, x, y, parms0, lower=[0.03, 0.04, 0.3, 2.35, -Inf, -Inf], upper=[0.06, 0.08, 0.35, 2.4, Inf, Inf])
params[n, :] = fit.param
Plots.plot(x_array, y)
Plots.plot!(x_array, [model(x, fit.param) for x in x_array])






for n in 1:partitions
    y = data_vec[n]
    C = sum(y) / length(y)
    A = 2 * sum(abs.(y .- C)) / length(y)
    B = 0
    parms0 = [2.0, A, B, C]

    x_array = Vector{Float64}(0:(6.28)/(length(y)-1):Float64(6.28))

    fitting = LsqFit.curve_fit(model, x_array, y, parms0; lower=[0.0, 0.25 * A, -π, 0.25 * C], upper=[Inf, 4 * A, Float64(π), 4 * C])

    params[n, :] = fitting.param
    # Plots.plot(x_array, y)
    # Plots.plot!(x_array, [model(x, fitting.param) for x in x_array])
end

@. linear_model(x, parms2) = parms2[2] + parms2[1] * x
p0 = [0.0, 0]
linear_fit = LsqFit.curve_fit(linear_model, Vector{Float64}(1:partitions), params[:, 3], p0)
slope = linear_fit.param[1]
error = LsqFit.estimate_errors(linear_fit)[1]









θ0 = (89) * (pi / 180)
range = 0.01
slope_vec = []
error_vec = []
θ_vec = Vector((θ0-range):0.001:(θ0+range))


for θ in θ_vec
    (N, M) = size(raw_image)
    LX = M - 1
    LY = N - 1
    #Aspect = 1.2
    #Aspect = 2.5
    Aspect = 10
    (R, S) = rotation_crop(dims=size(raw_image), angle=θ, AR=Aspect)


    rotated = OffsetArrays.centered(ImageTransformations.imrotate(raw_image, θ))



    sample = Float64.(rotated[-R:R, -S:S])



    strip_length = 10
    strip_skip = 1

    partitions = ((size(sample)[1] - strip_length) ÷ strip_skip) + 1
    partition_rngs = [(strip_skip*(n-1)+1):(strip_skip*(n-1)+strip_length) for n in 1:partitions]

    partition_rngs[length(partition_rngs)]
    data_vec = [dropdims(sum(sample[rng, :], dims=(1)), dims=(1)) / length(rng) for rng in partition_rngs]


    params = zeros(partitions, 4)
    model(r, parms::Vector{Float64}) = parms[2] * cos.(parms[1] * r .+ parms[3]) .+ parms[4]






    for n in 1:partitions
        y = data_vec[n]
        C = sum(y) / length(y)
        A = 2 * sum(abs.(y .- C)) / length(y)
        B = 0
        parms0 = [2.0, A, B, C]

        x_array = Vector{Float64}(0:(6.28)/(length(y)-1):Float64(6.28))

        fitting = LsqFit.curve_fit(model, x_array, y, parms0; lower=[0.0, 0.25 * A, -π, 0.25 * C], upper=[Inf, 4 * A, Float64(π), 4 * C])

        params[n, :] = fitting.param
        # Plots.plot(x_array, y)
        # Plots.plot!(x_array, [model(x, fitting.param) for x in x_array])
    end

    @. linear_model(x, parms2) = parms2[2] + parms2[1] * x
    p0 = [0.0, 0]
    linear_fit = LsqFit.curve_fit(linear_model, Vector{Float64}(1:partitions), params[:, 3], p0)
    slope = linear_fit.param[1]
    error = LsqFit.estimate_errors(linear_fit)[1]
    append!(slope_vec, slope)
    append!(error_vec, error)
end


error_vec

Plots.plot((180 / pi) * θ_vec, slope_vec, yerror=error_vec, legend=:none, xlabel="Angle (degrees)", ylabel="Phase Slope (Arbitrary Scale)", seriestype=:scatter)
#Plots.savefig("1f_aligned.pdf")


final_fit = LsqFit.curve_fit(linear_model, θ_vec, slope_vec, p0)
final_parms = final_fit.param
finalangle = -final_parms[2] / final_parms[1]
finalangle * (180 / pi)

