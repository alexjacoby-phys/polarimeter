import Images, FileIO, LinearAlgebra, DelimitedFiles
using Interpolations
function rescale(x::Float64; α::Float64=1.0)
    return (exp(α * x) - 1) / (exp(α) - 1)
end





ϕ_vec = Vector{Float64}([])




fn_base = string("//Volumes/AJACOBY/long_angle_benchmark/")










begin
    fn = string(fn_base, "20")
    raw_image = Images.Gray.(Images.load(string(fn, ".png")))[100:900, 300:1300]
    dat = Float64.(raw_image)
    (N, M) = size(dat)



    dat = dat - ones(Float64, size(dat)...) * sum(dat) / *(size(dat)...)





    N_X_MAX = 2.75
    N_X_MIN = 1.
    N_Y_MAX = 4
    N_Y_MIN = 2
    resolution = 0.001




    K_X = vcat(Vector{Float64}(2 * pi * (-N_X_MAX:resolution:-N_X_MIN) / M), Vector{Float64}(2 * pi * (N_X_MIN:resolution:N_X_MAX) / M))
    K_Y = vcat(Vector{Float64}(2 * pi * (-N_Y_MAX:resolution:-N_Y_MIN) / N), Vector{Float64}(2 * pi * (N_Y_MIN:resolution:N_Y_MAX) / N))


    pFTL = exp.(-im * [i * j for i in K_Y, j in 1:N])

    pFTR = exp.(-im * [i * j for i in 1:M, j in K_X])




    pftdat = abs.(pFTL * dat * pFTR)
    # A = max(pftdat...)
    # pftdat = pftdat * (1 / A)
    # downsampled = pftdat[1:10:size(pftdat)[1], 1:10:size(pftdat)[2]]
    # downsampled = downsampled / max(downsampled...)
    # Images.Gray.(downsampled)




    (ymax, xmax) = Tuple(findmax(pftdat)[2])
    (kx, ky) = (K_X[xmax], K_Y[ymax])


    ϕ = atan(ky, kx)
    ϕ = mod(ϕ + π / 2, π) - π / 2
    ϕ_deg = ϕ * (180 / pi)
    append!(ϕ_vec, ϕ_deg)
end


























begin
    fn = string(fn_base,"40")
    raw_image = Images.Gray.(Images.load(string(fn, ".png")))[100:900, 300:1300]
    dat = Float64.(raw_image)
    (N, M) = size(dat)



    dat = dat - ones(Float64, size(dat)...) * sum(dat) / *(size(dat)...)





    N_X_MAX = 4.0
    N_X_MIN = 2.1
    N_Y_MAX = 3.15
    N_Y_MIN = 1.3
    resolution = 0.001




    K_X = vcat(Vector{Float64}(2 * pi * (-N_X_MAX:resolution:-N_X_MIN) / M), Vector{Float64}(2 * pi * (N_X_MIN:resolution:N_X_MAX) / M))
    K_Y = vcat(Vector{Float64}(2 * pi * (-N_Y_MAX:resolution:-N_Y_MIN) / N), Vector{Float64}(2 * pi * (N_Y_MIN:resolution:N_Y_MAX) / N))


    pFTL = exp.(-im * [i * j for i in K_Y, j in 1:N])

    pFTR = exp.(-im * [i * j for i in 1:M, j in K_X])




    pftdat = abs.(pFTL * dat * pFTR)
    A = max(pftdat...)
    pftdat = pftdat*(1/A)
    downsampled = pftdat[1:10:size(pftdat)[1], 1:10:size(pftdat)[2]]
    downsampled = downsampled / max(downsampled...)
    Images.Gray.(downsampled)




    (ymax, xmax) = Tuple(findmax(pftdat)[2])
    (kx, ky) = (K_X[xmax], K_Y[ymax])


    ϕ = atan(ky, kx)
    ϕ = mod(ϕ + π / 2, π) - π / 2
    ϕ_deg = ϕ * (180 / pi)
    append!(ϕ_vec,ϕ_deg)
end






















begin
    fn = string(fn_base, "60")
    raw_image = Images.Gray.(Images.load(string(fn, ".png")))[100:900, 300:1300]
    dat = Float64.(raw_image)
    (N, M) = size(dat)



    dat = dat - ones(Float64, size(dat)...) * sum(dat) / *(size(dat)...)





    N_X_MAX = 4.75
    N_X_MIN = 2.75
    N_Y_MAX = 2.4
    N_Y_MIN = 0.4
    resolution = 0.001




    K_X = vcat(Vector{Float64}(2 * pi * (-N_X_MAX:resolution:-N_X_MIN) / M), Vector{Float64}(2 * pi * (N_X_MIN:resolution:N_X_MAX) / M))
    K_Y = vcat(Vector{Float64}(2 * pi * (-N_Y_MAX:resolution:-N_Y_MIN) / N), Vector{Float64}(2 * pi * (N_Y_MIN:resolution:N_Y_MAX) / N))


    pFTL = exp.(-im * [i * j for i in K_Y, j in 1:N])

    pFTR = exp.(-im * [i * j for i in 1:M, j in K_X])




    pftdat = abs.(pFTL * dat * pFTR)
    A = max(pftdat...)
    pftdat = pftdat * (1 / A)
    downsampled = pftdat[1:10:size(pftdat)[1], 1:10:size(pftdat)[2]]
    downsampled = downsampled / max(downsampled...)
    Images.Gray.(downsampled)




    (ymax, xmax) = Tuple(findmax(pftdat)[2])
    (kx, ky) = (K_X[xmax], K_Y[ymax])


    ϕ = atan(ky, kx)
    ϕ = mod(ϕ + π / 2, π) - π / 2
    ϕ_deg = ϕ * (180 / pi)
    append!(ϕ_vec, ϕ_deg)

end




































begin

    fn = string(fn_base, "80")
    raw_image = Images.Gray.(Images.load(string(fn, ".png")))[100:900, 300:1300]
    dat = Float64.(raw_image)
    (N, M) = size(dat)



    dat = dat - ones(Float64, size(dat)...) * sum(dat) / *(size(dat)...)





    N_X_MAX = 5.25
    N_X_MIN = 3.25
    N_Y_MAX = 1.25
    N_Y_MIN = 0.0
    resolution = 0.001




    K_X = vcat(Vector{Float64}(2 * pi * (-N_X_MAX:resolution:-N_X_MIN) / M), Vector{Float64}(2 * pi * (N_X_MIN:resolution:N_X_MAX) / M))
    K_Y = vcat(Vector{Float64}(2 * pi * (-N_Y_MAX:resolution:-N_Y_MIN) / N), Vector{Float64}(2 * pi * (N_Y_MIN:resolution:N_Y_MAX) / N))


    pFTL = exp.(-im * [i * j for i in K_Y, j in 1:N])

    pFTR = exp.(-im * [i * j for i in 1:M, j in K_X])




    pftdat = abs.(pFTL * dat * pFTR)
    A = max(pftdat...)
    pftdat = pftdat * (1 / A)
    downsampled = pftdat[1:10:size(pftdat)[1], 1:10:size(pftdat)[2]]
    downsampled = downsampled / max(downsampled...)
    Images.Gray.(downsampled)




    (ymax, xmax) = Tuple(findmax(pftdat)[2])
    (kx, ky) = (K_X[xmax], K_Y[ymax])


    ϕ = atan(ky, kx)
    ϕ = mod(ϕ + π / 2, π) - π / 2
    ϕ_deg = ϕ * (180 / pi)
    append!(ϕ_vec, ϕ_deg)

end







[ϕ_vec[i+1]-ϕ_vec[i] for i in 1:3]