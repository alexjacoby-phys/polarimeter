import Images, FileIO, LinearAlgebra, LsqFit, ImageTransformations, OffsetArrays, Meshes, Plots
include("/Users/alexjacoby/Documents/Research_Code/polarimeter/basic_functionality/Rotation_Rectangle.jl")
include("/Users/alexjacoby/Documents/Research_Code/polarimeter/basic_functionality/1d_averager.jl")


small_fn_vec = vcat(reverse([string("-", i) for i in 1:30]), "0", [string("+", i) for i in 1:20])
final_angle_vec = zeros(length(small_fn_vec))
final_angle_vec_deg = zeros(length(small_fn_vec))



for (α, small_fn) in pairs(small_fn_vec)
    γ = 12.6 * (π / 180)
    begin
        raw_imagepr = Images.Gray.(Images.load(string("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv4/peaks/peaks_png/",small_fn,".png")))

        raw_image = OffsetArrays.centered(ImageTransformations.imrotate(raw_imagepr, γ))[-240:280, 20:600]
        dat = Float64.(raw_image)

        correction_imagepr = Images.Gray.(Images.load("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv4/background_low_intensity.png"))
        correction_image = OffsetArrays.centered(ImageTransformations.imrotate(correction_imagepr, γ))[-240:280, 20:600]

        correction = Float64.(correction_image) .^ (-1)
        correction = (*(size(correction)...) / sum(correction)) * correction
        dat = dat .* correction
        image = Images.Gray.(dat)
    end








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
        Aspect = LX / LY
        (R, S) = rotation_crop(dims=size(image), angle=θ, AR=Aspect)
        rotated = OffsetArrays.centered(ImageTransformations.imrotate(image, θ))
        # Plots.plot(rotated);
        # Plots.plot!([S, -S, S, -S], [R, R, -R, -R], seriestype=:scatter, color=:red)


        sample = Float64.(rotated[-R:R, -S:S])



        strip_length = 200
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


        println(string("θ = ", θ, " done"))
    end




    Plots.plot((180 / pi) * (θ_vec .+ γ), slope_vec1, yerror=error_vec1, legend=:none, xlabel="Angle (degrees)", ylabel="Phase Slope (Arbitrary Scale)", seriestype=:scatter)
    Plots.savefig(string("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv4/peaks/Phase_Slope_Figs/", small_fn, ".pdf"))

    linear_model(x::Vector{Float64}, parms2::Vector{Float64}) = parms2[2] .+ parms2[1] * x
    p0 = [0.0, 0]
    final_fit = LsqFit.curve_fit(linear_model, θ_vec, slope_vec1, p0)
    final_parms = final_fit.param
    finalangle_unshifted = -final_parms[2] / final_parms[1]
    finalangle = (finalangle_unshifted + γ)
    finalangle_deg = (finalangle_unshifted + γ) * (180 / pi)
    final_angle_vec[α] = finalangle
    final_angle_vec_deg[α] = finalangle_deg

    # Plots.plot((180 / pi) * θ_vec, slope_vec1, yerror=error_vec1, legend=:none, xlabel="Angle (degrees)", ylabel="Phase Slope (Arbitrary Scale)", seriestype=:scatter);
    # Plots.plot!((180 / pi) * θ_vec, linear_model(θ_vec, final_parms))
end

DelimitedFiles.writedlm("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv4/peaks/exact_angles.txt", [small_fn_vec, final_angle_vec, final_angle_vec_deg])

Plots.plot(final_angle_vec_deg)