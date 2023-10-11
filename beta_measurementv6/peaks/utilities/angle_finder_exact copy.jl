import Images, FileIO, LinearAlgebra, LsqFit, ImageTransformations, OffsetArrays, Meshes, Plots
include("/Users/alexjacoby/Documents/Research_Code/polarimeter/basic_functionality/Rotation_Rectangle.jl")
include("/Users/alexjacoby/Documents/Research_Code/polarimeter/basic_functionality/1d_averager.jl")



small_fn_vec = small_fn_vec = vcat(reverse([string("-", i) for i in 1:25]), "0", [string("+", i) for i in 1:25])


n = 26

small_fn = small_fn_vec[n]
γ = -(66.5) * (π / 180)
# begin
    raw_imagepr = Images.Gray.(Images.load(string("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv6/peaks/peaks_png/",small_fn,".png")))

    raw_image = OffsetArrays.centered(ImageTransformations.imrotate(raw_imagepr, γ))[-450:0, -250:250]
    dat = Float64.(raw_image)

    correction_imagepr = Images.Gray.(Images.load("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv6/background/backgroundv3.png"))
    correction_image = OffsetArrays.centered(ImageTransformations.imrotate(correction_imagepr, γ))[-450:0, -250:250]

    correction = Float64.(correction_image) .^(-1)
    correction = (*(size(correction)...) / sum(correction)) * correction;
    dat = dat .* correction
    image = Images.Gray.(dat)
# end




# using Images, FFTW, LinearAlgebra

# function fermi(E::Number; β::Number, μ::Number)
#     return 1 / (exp(β * (E - μ)) + 1)
# end




# FT_DAT = FFTW.fftshift(FFTW.fft(dat))


# (N, M) = size(dat)
# distances = [LinearAlgebra.norm([y - N / 2, x - M / 2]) for y in 1:N, x in 1:M]
# distances = distances / max(distances...)

# ffilter = fermi.(distances, β=20.0, μ=0.05)

# Images.Gray.(abs.(FT_DAT) / (sqrt(N) * sqrt(M)))
# Images.Gray.(ffilter)


# Images.Gray.(abs.(FFTW.ifft(FFTW.fftshift(FT_DAT .* ffilter))))

# dat = abs.(FFTW.ifft(FFTW.fftshift(FT_DAT .* ffilter)))








θ = 0. *(π/180) #1.48 *(π / 180)
(N, M) = size(image)
LX = M - 1
LY = N - 1
#Aspect = 1.2
#Aspect = 2.5
Aspect = LX/LY
(R, S) = rotation_crop(dims=size(image), angle=θ, AR=Aspect)
rotated = OffsetArrays.centered(ImageTransformations.imrotate(image, θ));
# Plots.plot(rotated);
# Plots.plot!([S, -S, S, -S], [R, R, -R, -R], seriestype=:scatter, color=:red)


sample = Float64.(rotated[-R:R, -S:S])
Images.Gray.(sample);


strip_length = 40
strip_skip = 2

partitions = ((size(sample)[1] - strip_length) ÷ strip_skip) + 1
partition_rngs = [(strip_skip*(n-1)+1):(strip_skip*(n-1)+strip_length) for n in 1:partitions]

partition_rngs[length(partition_rngs)]
data_vec = [dropdims(sum(sample[rng, :], dims=(1)), dims=(1)) / length(rng) for rng in partition_rngs]





model(r::Vector{Float64}, parms::Vector{Float64}) = parms[1] * cos.(parms[4] * r .- parms[5]) .+ parms[2] * cos.(2 * parms[4] * r .+ parms[6] .- 2*parms[5]) .+ parms[3]



params = zeros(partitions,6)


parms0 = [0.1, 0.05, 0.2, 2.5, 0, 0]
for n in 1:partitions
    y = data_vec[n]
    x = π * Vector(-1.0:2/(length(y)-1):1)
    parms0 = [0.1, 0.05, 0.2, 2.5, 0, 0]

    fit = LsqFit.curve_fit(model, x, y, parms0)
    params[n, :] = cf_2f(fit.param)
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
Plots.plot!(1:partitions, linear_model(Vector(1:partitions),linear_fit.param))




animation = Plots.@animate for n in 1:partitions
    parms0 = params[n,:]
    y = data_vec[n]
    x = π * Vector(-1.0:2/(length(y)-1):1)
    fit = LsqFit.curve_fit(model, x, y, parms0)
    params[n, :] = fit.param
    Plots.plot(x, y, label="Data", seriestype=:scatter, legend=:topright)
    Plots.plot!(x, model(x, fit.param), label="Fitted", linewidth=5)
end




Plots.gif(animation, "moving.gif", fps=15)













θ0 = 0 * (pi / 180)
range = 2 * (pi / 180)
increments = 20
slope_vec1 = []
error_vec1 = []
θ_vec = Vector((θ0-range):(2*range/(increments-1)):(θ0+range))





for θ in θ_vec
    (N, M) = size(image)
    LX = M - 1
    LY = N - 1
    #Aspect = 1.2
    #Aspect = 2.5
    Aspect = LX/LY
    (R, S) = rotation_crop(dims=size(image), angle=θ, AR=Aspect)
    rotated = OffsetArrays.centered(ImageTransformations.imrotate(image, θ))
    # Plots.plot(rotated);
    # Plots.plot!([S, -S, S, -S], [R, R, -R, -R], seriestype=:scatter, color=:red)


    sample = Float64.(rotated[-R:R, -S:S])



    strip_length = 40
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


    linear_fit2 = LsqFit.curve_fit(linear_model, Vector{Float64}(1:partitions), params[:, 6], p0)
    slope2 = linear_fit2.param[1]
    error2 = LsqFit.estimate_errors(linear_fit2)[1]


    append!(slope_vec1, slope1)
    append!(error_vec1, error1)


    println(string("θ = ",θ," done"))
end


slope_vec1
θ_vec


Plots.plot((180 / pi) * (θ_vec .+ γ), slope_vec1, yerror=error_vec1, legend=:none, xlabel="Angle (degrees)", ylabel="Phase Slope (Arbitrary Scale)", seriestype=:scatter)

θ_vec
slope_vec1
linear_model(x::Vector{Float64}, parms2::Vector{Float64}) = parms2[2] .+ parms2[1] * x
p0 = [0.0, 0]
final_fit = LsqFit.curve_fit(linear_model, θ_vec, slope_vec1, p0)
final_parms = final_fit.param
finalangle = -final_parms[2] / final_parms[1]
(finalangle + γ) * (180 / pi)

Plots.plot((180 / pi) * (θ_vec .+ γ), slope_vec1, yerror=error_vec1, legend=:none, xlabel="Angle (degrees)", ylabel="Phase Slope (Arbitrary Scale)", seriestype=:scatter);
Plots.plot!((180 / pi) * (θ_vec .+γ), linear_model( θ_vec,final_parms))





