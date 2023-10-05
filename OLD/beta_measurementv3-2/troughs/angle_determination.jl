import Images, FileIO, LinearAlgebra, LsqFit, ImageTransformations, OffsetArrays, Meshes



raw_image = Images.Gray.(Images.load("//Volumes/AJACOBY/beta_measurementv3/troughs/-7.png"))[100:900, 300:1300]

dat = Float64.(raw_image)

correction = (Float64.(Images.Gray.(Images.load("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv3-2/intensity_correction.png"))[100:900, 300:1300])) .^ (-1)
correction = (*(size(correction)...) / sum(correction)) * correction;
dat = dat .* correction
Images.Gray.(dat)
function rotation(angle::Float64)
    return [cos(angle) -sin(angle); sin(angle) cos(angle)]
end



θ = -27.177 * (pi / 180)


(N, M) = size(raw_image)
LX = M - 1
LY = N - 1






Aspect = 1.
(R, S) = rotation_crop(dims=size(raw_image), angle=θ, AR= Aspect)


rotated = OffsetArrays.centered(ImageTransformations.imrotate(raw_image, θ))

#rotated = OffsetArrays.centered(ImageTransformations.imrotate(raw_image, -final_parms[2] / final_parms[1]))
Plots.plot(rotated);
Plots.plot!([S, -S, S, -S], [R, R, -R, -R], seriestype=:scatter, color=:red)


sample = Float64.(rotated[-R:R, -S:S])
Images.Gray.(sample)


strip_length = 20
strip_skip = 1

partitions = ((size(sample)[1] - strip_length) ÷ strip_skip) + 1
partition_rngs = [(strip_skip*(n-1)+1):(strip_skip*(n-1)+strip_length) for n in 1:partitions]

partition_rngs[length(partition_rngs)]
data_vec = [dropdims(sum(sample[rng, :], dims=(1)), dims=(1)) / length(rng) for rng in partition_rngs]


params = zeros(partitions, 6)
model(r, parms::Vector{Float64}) = parms[2] * cos.(parms[1] * r .+ parms[4]) .+ parms[3] * cos.(2 * parms[1] * r .+ parms[5]) .+ parms[6]

n = 1
y = data_vec[n]
C = sum(y) / length(y)
A = 2 * sum(abs.(y .- C)) / length(y)

parms0 = [4.0, A, A, 0, 0, C]

x_array = Vector{Float64}(0:(6.28)/(length(y)-1):Float64(6.28))

fitting = LsqFit.curve_fit(model, x_array, y, parms0; lower=[0.0, 0.0, 0.0, -π, -π, 0.0], upper=[Inf, Inf, Inf, π, π, Inf])

params[n, :] = fitting.param
Plots.plot(x_array, y)
Plots.plot!(x_array, [model(x, fitting.param) for x in x_array])




for n in 1:partitions
    y = data_vec[n]
    C = sum(y) / length(y)
    A = 2 * sum(abs.(y .- C)) / length(y)

    parms0 = [4.0, A, A, 0, 0, C]

    x_array = Vector{Float64}(0:(6.28)/(length(y)-1):Float64(6.28))

    fitting = LsqFit.curve_fit(model, x_array, y, parms0; lower=[0.0, 0.0, 0.0, -π, -π, 0.0], upper=[Inf, Inf, Inf, π, π, Inf])

    params[n, :] = fitting.param
end

@. linear_model(x, parms2) = parms2[2] + parms2[1] * x
p0 = [0.0, 0]
linear_fit = LsqFit.curve_fit(linear_model, Vector{Float64}(1:partitions), params[:, 3], p0)
slope = linear_fit.param[1]



Plots.plot(1:partitions, params[:, 3])
Plots.plot!(1:partitions, [linear_model(k, linear_fit.param) for k in 1:partitions])











θ0 = -23* (pi / 180)
range = 0.01
slope_vec = []
error_vec = []
θ_vec = Vector((θ0-range):0.001:(θ0+range))
for θ in θ_vec
    #θ = -0.152
    (N, M) = size(raw_image)
    LX = M - 1
    LY = N - 1


    (PX, PY) = ((LX * cos(θ) - LY * sin(θ)) / 2, (LY * cos(θ) + LX * sin(θ)) / 2)



    AR = 2.5
    (R, S) = rotation_crop(dims=size(raw_image), angle=θ, AR=2.5)


    rotated = OffsetArrays.centered(ImageTransformations.imrotate(raw_image, θ))



    sample = Float64.(rotated[-R:R, -S:S])



    strip_length = 10
    strip_skip = 1

    partitions = ((size(sample)[1] - strip_length) ÷ strip_skip) + 1
    partition_rngs = [(strip_skip*(n-1)+1):(strip_skip*(n-1)+strip_length) for n in 1:partitions]

    partition_rngs[length(partition_rngs)]
    data_vec = [dropdims(sum(sample[rng, :], dims=(1)), dims=(1)) / length(rng) for rng in partition_rngs]


    params = zeros(partitions, 6)
    model(r, parms::Vector{Float64}) = parms[2] * cos.(parms[1] * r .+ parms[4]) .+ parms[3] * cos.(2 * parms[1] * r .+ parms[5]) .+ parms[6]






    for n in 1:partitions
        y = data_vec[n]
        C = sum(y) / length(y)
        A = 2 * sum(abs.(y .- C)) / length(y)

        parms0 = [4.0, A, A, 0, 0, C]

        x_array = Vector{Float64}(0:(6.28)/(length(y)-1):Float64(6.28))

        fitting = LsqFit.curve_fit(model, x_array, y, parms0; lower=[0.0, 0.0, 0.0, -π, -π, 0.0], upper=[Inf, Inf, Inf, π, π, Inf])

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

Plots.plot((180 / pi) * θ_vec, slope_vec, yerror=error_vec, legend=:none, xlabel="Angle (degrees)", ylabel="Phase Slope (Arbitrary Scale)", seriestype=:scatter)



final_fit = LsqFit.curve_fit(linear_model, θ_vec, slope_vec, p0)
final_parms = final_fit.param
-final_parms[2] / final_parms[1]






































































raw_image = Images.Gray.(100 * Images.load("/Users/alexjacoby/Documents/Research_Code/polarimeter/1faligned.tif"))[270:710, 700:1310]

dat = Float64.(raw_image)

correction = (Float64.(Images.Gray.(Images.load("/Users/alexjacoby/Documents/Research_Code/polarimeter/backgroundnowedge.tif"))[270:710, 700:1310])) .^ (-1)
correction = (*(size(correction)...) / sum(correction)) * correction;
dat = dat .* correction
Images.Gray.(dat)
function rotation(angle::Float64)
    return [cos(angle) -sin(angle); sin(angle) cos(angle)]
end



θ = 88.0 * (pi / 180)#Float64(π/2 - 0.03)




(N, M) = size(raw_image)
LX = M - 1
LY = N - 1


Aspect = 1
(R, S) = rotation_crop(dims=size(raw_image), angle=θ, AR=Aspect)


rotated = OffsetArrays.centered(ImageTransformations.imrotate(raw_image, θ));



sample = Float64.(rotated[-R:R, -S:S])
Images.Gray.(sample)


strip_length = 10
strip_skip = 1

partitions = ((size(sample)[1] - strip_length) ÷ strip_skip) + 1
partition_rngs = [(strip_skip*(n-1)+1):(strip_skip*(n-1)+strip_length) for n in 1:partitions]

partition_rngs[length(partition_rngs)]
data_vec = [dropdims(sum(sample[rng, :], dims=(1)), dims=(1)) / length(rng) for rng in partition_rngs]


params = zeros(partitions, 4)
model(r, parms::Vector{Float64}) = parms[2] * cos.(parms[1] * r .+ parms[3]) .+ parms[4]




n = 400
y = data_vec[n]
C = sum(y) / length(y)
A = 2 * sum(abs.(y .- C)) / length(y)
B = 0
parms0 = [2.0, A, B, C]

x_array = Vector{Float64}(0:(6.28)/(length(y)-1):Float64(6.28))

fitting = LsqFit.curve_fit(model, x_array, y, parms0; lower=[0.0, 0.25 * A, -π, 0.25 * C], upper=[Inf, 4 * A, Float64(π), 4 * C])

params[n, :] = fitting.param
Plots.plot(x_array, y);
Plots.plot!(x_array, [model(x, fitting.param) for x in x_array])





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

Plots.plot(1:partitions, params[:, 3])

















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

































