import Images, FileIO, LinearAlgebra, LsqFit, ImageTransformations, OffsetArrays, Meshes
include("/Users/alexjacoby/Documents/Research_Code/polarimeter/basic_functionality/Rotation_Rectangle.jl")
include("/Users/alexjacoby/Documents/Research_Code/polarimeter/basic_functionality/1d_averager.jl")



raw_image = Images.Gray.(Images.load("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv4/1f/1f_png/+2.png"))[270:710, 700:1350]

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

filter = fermi.(distances, β=20.0, μ=0.05)

Images.Gray.(abs.(FT_DAT) / (sqrt(N) * sqrt(M)))
Images.Gray.(filter)


Images.Gray.(abs.(FFTW.ifft(FFTW.fftshift(FT_DAT .* filter))))

dat = abs.(FFTW.ifft(FFTW.fftshift(FT_DAT .* filter)))





θ = -9. * (π / 180)
(N, M) = size(raw_image)
LX = M - 1
LY = N - 1
#Aspect = 1.2
#Aspect = 2.5
Aspect = 1.5
(R, S) = rotation_crop(dims=size(raw_image), angle=θ, AR=Aspect)


rotated = OffsetArrays.centered(ImageTransformations.imrotate(raw_image, θ));
Plots.plot(rotated);
Plots.plot!([S, -S, S, -S], [R, R, -R, -R], seriestype=:scatter, color=:red)


sample = Float64.(rotated[-R:R, -S:S])



strip_length = 100
strip_skip = 1

partitions = ((size(sample)[1] - strip_length) ÷ strip_skip) + 1
partition_rngs = [(strip_skip*(n-1)+1):(strip_skip*(n-1)+strip_length) for n in 1:partitions]

partition_rngs[length(partition_rngs)]
data_vec = [dropdims(sum(sample[rng, :], dims=(1)), dims=(1)) / length(rng) for rng in partition_rngs]


model(r::Vector{Float64}, parms::Vector{Float64}) = parms[1] * cos.(parms[3] * r .+ parms[4]) .+ parms[2]
params = zeros(partitions, 4)



y = data_vec[n]
x = π * Vector(-1.0:2/(length(y)-1):1)
parms0 = [0.25, 0.25, 2.0, π]

fit = LsqFit.curve_fit(model, x, y, parms0)#, lower=[0.0, 0.0, 0.0, 1., -π, -π], upper=[0.5, 0.5, 0.5, 3.5, π, π])
params[n, :] = fit.param
Plots.plot(x, y);

Plots.plot!(x, model(x, fit.param))


for n in 1:partitions
    y = data_vec[n]
    x = π * Vector(-1.0:2/(length(y)-1):1)
    parms0 = [0.25, 0.25, 2.0, π]

    fit = LsqFit.curve_fit(model, x, y, parms0)#, lower=[0.0, 0.0, 0.0, 1., -π, -π], upper=[0.5, 0.5, 0.5, 3.5, π, π])
    params[n, :] = fit.param
    Plots.plot(x, y)
    Plots.plot!(x, model(x, parms0))
    Plots.plot!(x, model(x, fit.param))
end

Plots.plot(1:partitions, mod.(params[:, 4], 2π))



linear_model(x::Vector{Float64}, parms2::Vector{Float64}) = parms2[2] .+ parms2[1] * x
p0 = [0.0, 0]
linear_fit = LsqFit.curve_fit(linear_model, Vector{Float64}(1:partitions), params[:, 4], p0)
slope = linear_fit.param[1]
error = LsqFit.estimate_errors(linear_fit)[1]


Plots.plot(1:partitions, mod.(params[:, 4], 2π));
Plots.plot!(1:partitions, linear_model(1:partitions, linear_fit.param))




















θ0 = -9. *(pi/180)
range = 0.01
slope_vec = []
error_vec = []
θ_vec = Vector((θ0-range):0.001:(θ0+range))
for θ in θ_vec
    Aspect = 1.5
    (R, S) = rotation_crop(dims=size(raw_image), angle=θ, AR=Aspect)


    rotated = OffsetArrays.centered(ImageTransformations.imrotate(raw_image, θ))
    # Plots.plot(rotated)
    # Plots.plot!([S, -S, S, -S], [R, R, -R, -R], seriestype=:scatter, color=:red)


    sample = Float64.(rotated[-R:R, -S:S])



    strip_length = 100
    strip_skip = 1

    partitions = ((size(sample)[1] - strip_length) ÷ strip_skip) + 1
    partition_rngs = [(strip_skip*(n-1)+1):(strip_skip*(n-1)+strip_length) for n in 1:partitions]

    partition_rngs[length(partition_rngs)]
    data_vec = [dropdims(sum(sample[rng, :], dims=(1)), dims=(1)) / length(rng) for rng in partition_rngs]


    model(r::Vector{Float64}, parms::Vector{Float64}) = parms[1] * cos.(parms[3] * r .+ parms[4]) .+ parms[2]
    params = zeros(partitions, 4)



    for n in 1:partitions
        y = data_vec[n]
        x = π * Vector(-1.0:2/(length(y)-1):1)
        parms0 = [0.25, 0.25, 2.0, π]

        fit = LsqFit.curve_fit(model, x, y, parms0)#, lower=[0.0, 0.0, 0.0, 1., -π, -π], upper=[0.5, 0.5, 0.5, 3.5, π, π])
        params[n, :] = fit.param
        # Plots.plot(x, y)
        # Plots.plot!(x, model(x, parms0))
        # Plots.plot!(x, model(x, fit.param))
    end


    linear_model(x::Vector{Float64}, parms2::Vector{Float64}) = parms2[2] .+ parms2[1] * x
    p0 = [0.0, 0]
    linear_fit = LsqFit.curve_fit(linear_model, Vector{Float64}(1:partitions), params[:, 4], p0)
    slope = linear_fit.param[1]
    error = LsqFit.estimate_errors(linear_fit)[1]

    append!(slope_vec, slope)
    append!(error_vec, error)
end

Plots.plot((180 / pi) * θ_vec, slope_vec, yerror=error_vec, legend=:none, xlabel="Angle (degrees)", ylabel="Phase Slope (Arbitrary Scale)", seriestype=:scatter)



final_fit = LsqFit.curve_fit(linear_model, θ_vec, slope_vec, p0)
final_parms = final_fit.param
(-final_parms[2] / final_parms[1])*180/pi



index = DelimitedFiles.readdlm("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv4/1f/Estimated_Angles.txt")
index[3,23]