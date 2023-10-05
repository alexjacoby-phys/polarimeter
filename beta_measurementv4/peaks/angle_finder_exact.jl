import Images, FileIO, LinearAlgebra, LsqFit, ImageTransformations, OffsetArrays, Meshes
include("/Users/alexjacoby/Documents/Research_Code/polarimeter/basic_functionality/Rotation_Rectangle.jl")
include("/Users/alexjacoby/Documents/Research_Code/polarimeter/basic_functionality/1d_averager.jl")



raw_image = Images.Gray.(Images.load("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv4/peaks/peaks_png/-10.png"))[270:710, 700:1350];

dat = Float64.(raw_image)

correction = (Float64.(Images.Gray.(Images.load("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv4/background_low_intensity.png"))[270:710, 700:1350])) .^ (-1)
correction = (*(size(correction)...) / sum(correction)) * correction;
dat = dat .* correction
Images.Gray.(dat)
function rotation(angle::Float64)
    return [cos(angle) -sin(angle); sin(angle) cos(angle)]
end




using Images, FFTW, LinearAlgebra

function fermi(E::Number; β::Number, μ::Number)
    return 1 / (exp(β * (E - μ)) + 1)
end




FT_DAT = FFTW.fftshift(FFTW.fft(dat))


(N, M) = size(dat)
distances = [LinearAlgebra.norm([y - N / 2, x - M / 2]) for y in 1:N, x in 1:M]
distances = distances / max(distances...)

ffilter = fermi.(distances, β=20.0, μ=0.05)

Images.Gray.(abs.(FT_DAT) / (sqrt(N) * sqrt(M)))
Images.Gray.(filter)


Images.Gray.(abs.(FFTW.ifft(FFTW.fftshift(FT_DAT .* ffilter))))

dat = abs.(FFTW.ifft(FFTW.fftshift(FT_DAT .* ffilter)))








θ = 14.3 * (π / 180)
(N, M) = size(raw_image)
LX = M - 1
LY = N - 1
#Aspect = 1.2
#Aspect = 2.5
Aspect =1.5
(R, S) = rotation_crop(dims=size(raw_image), angle=θ, AR=Aspect)
rotated = OffsetArrays.centered(ImageTransformations.imrotate(raw_image, θ))
Plots.plot(rotated);
Plots.plot!([S, -S, S, -S], [R, R, -R, -R], seriestype=:scatter, color=:red)


sample = Float64.(rotated[-R:R, -S:S])



strip_length = 100
strip_skip = 1

partitions = ((size(sample)[1] - strip_length) ÷ strip_skip) + 1
partition_rngs = [(strip_skip*(n-1)+1):(strip_skip*(n-1)+strip_length) for n in 1:partitions]

partition_rngs[length(partition_rngs)]
data_vec = [dropdims(sum(sample[rng, :], dims=(1)), dims=(1)) / length(rng) for rng in partition_rngs]





model(r::Vector{Float64}, parms::Vector{Float64}) = parms[1] * cos.(parms[4] * r .+ parms[5]) .+ parms[2] * cos.(2 * parms[4] * r .+ parms[6]) .+ parms[3]



params = zeros(partitions,6)


n=100
y = data_vec[n]
x = π * Vector(-1.0:2/(length(y)-1):1)
parms0 = [0.1, 0.05, 0.2, 2.0, 0, 0]

fit = LsqFit.curve_fit(model, x, y, parms0)#, lower=[0.0, 0.0, 0.0, 1., -π, -π], upper=[0.5, 0.5, 0.5, 3.5, π, π])
params[n, :] = fit.param
Plots.plot(x, y);

Plots.plot!(x, model(x, fit.param))



for n in 1:partitions
    y = data_vec[n]
    x = π * Vector(-1.0:2/(length(y)-1):1)
    parms0 = [0.1, 0.05, 0.2, 2.0, 0, 0]

    fit = LsqFit.curve_fit(model, x, y, parms0)#, lower=[0.0, 0.0, 0.0, 0., -Inf, -Inf], upper=[0.5, 0.5, 0.5, Inf, Inf, Inf])
    params[n, :] = fit.param
    # Plots.plot(x, y)
    # Plots.plot!(x, model(x, fit.param))
end

Plots.plot(1:partitions, mod.(params[:,5],π))
Plots.plot(1:partitions, mod.(params[:,6],2π))


@. linear_model(x, parms2) = parms2[2] + parms2[1] * x
p0 = [0.0, 0]
linear_fit = LsqFit.curve_fit(linear_model, Vector{Float64}(1:partitions), params[:, 3], p0)
slope = linear_fit.param[1]
error = LsqFit.estimate_errors(linear_fit)[1]









θ0 = (14.32) * (pi / 180)
range = 0.01
slope_vec1 = []
error_vec1 = []
slope_vec2 = []
error_vec2 = []
θ_vec = Vector((θ0-range):0.001:(θ0+range))


for θ in θ_vec
    (N, M) = size(raw_image)
    LX = M - 1
    LY = N - 1
    #Aspect = 1.2
    #Aspect = 2.5
    Aspect = 1.5
    (R, S) = rotation_crop(dims=size(raw_image), angle=θ, AR=Aspect)
    rotated = OffsetArrays.centered(ImageTransformations.imrotate(raw_image, θ))
    # Plots.plot(rotated);
    # Plots.plot!([S, -S, S, -S], [R, R, -R, -R], seriestype=:scatter, color=:red)


    sample = Float64.(rotated[-R:R, -S:S])



    strip_length = 100
    strip_skip = 1

    partitions = ((size(sample)[1] - strip_length) ÷ strip_skip) + 1
    partition_rngs = [(strip_skip*(n-1)+1):(strip_skip*(n-1)+strip_length) for n in 1:partitions]

    partition_rngs[length(partition_rngs)]
    data_vec = [dropdims(sum(sample[rng, :], dims=(1)), dims=(1)) / length(rng) for rng in partition_rngs]





    model(r::Vector{Float64}, parms::Vector{Float64}) = parms[1] * cos.(parms[4] * r .+ parms[5]) .+ parms[2] * cos.(2 * parms[4] * r .+ parms[6]) .+ parms[3]



    params = zeros(partitions, 6)


    for n in 1:partitions
        for n in 1:partitions
            y = data_vec[n]
            x = π * Vector(-1.0:2/(length(y)-1):1)
            parms0 = [0.1, 0.05, 0.2, 2.0, 0, 0]

            fit = LsqFit.curve_fit(model, x, y, parms0)#, lower=[0.0, 0.0, 0.0, 0., -Inf, -Inf], upper=[0.5, 0.5, 0.5, Inf, Inf, Inf])
            params[n, :] = fit.param
            # Plots.plot(x, y)
            # Plots.plot!(x, model(x, fit.param))
        end
    end



    linear_model(x::Vector{Float64}, parms2::Vector{Float64}) = parms2[2] .+ parms2[1] * x
    p0 = [0.0, 0]
    linear_fit1 = LsqFit.curve_fit(linear_model, Vector{Float64}(1:partitions), params[:, 5], p0)
    slope1 = linear_fit1.param[1]
    error1 = LsqFit.estimate_errors(linear_fit1)[1]


    linear_fit2 = LsqFit.curve_fit(linear_model, Vector{Float64}(1:partitions), params[:, 6], p0)
    slope2 = linear_fit2.param[1]
    error2 = LsqFit.estimate_errors(linear_fit2)[1]


    append!(slope_vec1, slope1)
    append!(error_vec1, error1)
    append!(slope_vec2, slope2)
    append!(error_vec2, error2)

    println(string("θ = ",θ," done"))
end




Plots.plot((180 / pi) * θ_vec, slope_vec, yerror=error_vec, legend=:none, xlabel="Angle (degrees)", ylabel="Phase Slope (Arbitrary Scale)", seriestype=:scatter)



final_fit = LsqFit.curve_fit(linear_model, θ_vec, slope_vec, p0)
final_parms = final_fit.param
finalangle = -final_parms[2] / final_parms[1]
finalangle * (180 / pi)










θ_vec = Vector(16:0.01:20)*π/180
σ = zeros(length(θ_vec))

for (n,θ) in pairs(θ_vec)
    (x,y) = convert_1d(dat, ϕ =θ)
    x = x[250:1750]
    y = y[250:1750]

    μ = sum(y)/length(y)
    σ[n]  = sum((y .- μ) .^2)
    println(θ)
end




Plots.plot(16:0.01:20,σ)