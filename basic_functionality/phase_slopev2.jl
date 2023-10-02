import Images, FileIO, LinearAlgebra, LsqFit, ImageTransformations, OffsetArrays



raw_image = Images.Gray.(Images.load("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv4/1f/-1.png"))[270:710, 700:1350];

dat = Float64.(raw_image)


correction = (Float64.(Images.Gray.(Images.load("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv4/FFiltered_Background_Low.png"))[270:710, 700:1350])) .^ (-1)
correction = (*(size(correction)...) / sum(correction)) * correction;
dat = dat .* correction
Images.Gray.(dat)
function rotation(angle::Float64)
    return [cos(angle) -sin(angle); sin(angle) cos(angle)]
end

θ0 = -0.16
range = 0.01
slope_vec = []

θ_vec = Vector((θ0-range):0.001:(θ0+range))
for θ in θ_vec
    #θ = -0.152
    (N, M) = size(raw_image)
    LX = M - 1
    LY = N - 1


    (PX, PY) = ((LX * cos(θ) - LY * sin(θ)) / 2, (LY * cos(θ) + LX * sin(θ)) / 2)



    AR = 0.6
    ψ = acos(1 / sqrt(1 + AR^2))
    xval = (PX * tan(-θ) + PY) / (tan(ψ) + abs(tan(θ)))
    R = Int64(floor((PX * tan(-θ) + PY) / (tan(ψ) + abs(tan(θ)))))
    S = Int64(floor(xval * AR))

    rotated = OffsetArrays.centered(ImageTransformations.imrotate(raw_image, θ))



    sample = Float64.(rotated[-S:S, -R:R])



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
        parms0 = [2.5, A, B, C]

        x_array = Vector{Float64}(0:(6.28)/(length(y)-1):Float64(6.28))

        fitting = LsqFit.curve_fit(model, x_array, y, parms0; lower=[0.0, 0.0, -π, 0.0], upper=[Inf, Inf, Float64(π), Inf])

        params[n, :] = fitting.param
        Plots.plot(x_array, y)
        Plots.plot!(x_array, [model(x, fitting.param) for x in x_array])
    end

    @. linear_model(x, parms2) = parms2[2] + parms2[1] * x
    p0 = [0.0, 0]
    linear_fit = LsqFit.curve_fit(linear_model, Vector{Float64}(1:partitions), params[:, 3], p0)
    slope = linear_fit.param[1]
    append!(slope_vec, slope)
end

Plots.plot(θ_vec, slope_vec)

final_fit = LsqFit.curve_fit(linear_model, θ_vec,slope_vec, p0)
final_parms = final_fit.param
-final_parms[2]/final_parms[1]












