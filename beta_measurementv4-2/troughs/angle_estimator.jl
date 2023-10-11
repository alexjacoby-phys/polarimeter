
using LinearAlgebra, Images, FileIO, DelimitedFiles











small_fn_vec = vcat(reverse([string("-", i) for i in 1:20]), "0", [string("+", i) for i in 1:20])
ϕ_vec = zeros(length(small_fn_vec))
ϕ_vec_deg = zeros(length(small_fn_vec))

for (ind, small_fn) in pairs(small_fn_vec)
    cd("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv4/troughs/troughs_png/")


    raw_image = Images.Gray.(Images.load(string(small_fn, ".png")))[270:710, 700:1350]
    dat = Float64.(raw_image)

    correction = (Float64.(Images.Gray.(Images.load("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv4/background_low_intensity.png"))[270:710, 700:1350])) .^ (-1)
    correction = (*(size(correction)...) / sum(correction)) * correction
    dat = dat .* correction








    (N, M) = size(dat)

    #dat = [0.5*(1+cos((i+j)/100)) for i in 1:1080, j in 1:1440]


    dat = dat - ones(Float64, size(dat)...) * sum(dat) / *(size(dat)...)




    N_X_MAX = 5.5
    N_X_MIN = 3
    N_Y_MAX = 3
    N_Y_MIN = 1
    resolution = 0.001




    K_X = vcat(Vector{Float64}(2 * pi * (-N_X_MAX:resolution:-N_X_MIN) / M), Vector{Float64}(2 * pi * (N_X_MIN:resolution:N_X_MAX) / M))
    K_Y = vcat(Vector{Float64}(2 * pi * (-N_Y_MAX:resolution:-N_Y_MIN) / N), Vector{Float64}(2 * pi * (N_Y_MIN:resolution:N_Y_MAX) / N))


    pFTL = exp.(-im * [i * j for i in K_Y, j in 1:N])

    pFTR = exp.(-im * [i * j for i in 1:M, j in K_X])







    pftdat = abs.(pFTL * dat * pFTR)
    # A = max(pftdat...)
    # pftdat = pftdat*(1/A)
    # downsampled = pftdat[1:10:size(pftdat)[1], 1:10:size(pftdat)[2]]
    # downsampled = downsampled / max(downsampled...)
    # Images.Gray.(downsampled)




    (ymax, xmax) = Tuple(findmax(pftdat)[2])
    (kx, ky) = (K_X[xmax], K_Y[ymax])


    ϕ = atan(ky, kx)
    ϕ = mod(ϕ + π / 2, π) - π / 2
    ϕ_deg = ϕ * (180 / pi)
    ϕ_vec[ind] = ϕ
    ϕ_vec_deg[ind] = ϕ_deg
end


cd("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv4/troughs")
DelimitedFiles.writedlm("Estimated_Angles.txt", [small_fn_vec, ϕ_vec, ϕ_vec_deg])


Plots.plot(ϕ_vec_deg)





# small_fn_vec = vcat(reverse([string("-", i) for i in 1:20]), "0", [string("+", i) for i in 1:20])
# fn = small_fn_vec[20]

# cd("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv4/troughs/troughs_png/")


# raw_image = Images.Gray.(Images.load(string(small_fn, ".png")))[270:710, 700:1350]
# dat = Float64.(raw_image)

# correction = (Float64.(Images.Gray.(Images.load("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv4/background_low_intensity.png"))[270:710, 700:1350])) .^ (-1)
# correction = (*(size(correction)...) / sum(correction)) * correction;
# dat = dat .* correction








# (N, M) = size(dat)

# #dat = [0.5*(1+cos((i+j)/100)) for i in 1:1080, j in 1:1440]


# dat = dat - ones(Float64, size(dat)...) * sum(dat) / *(size(dat)...)




# N_X_MAX = 5.5
# N_X_MIN = 3
# N_Y_MAX = 3
# N_Y_MIN = 1
# resolution = 0.01




# K_X = vcat(Vector{Float64}(2 * pi * (-N_X_MAX:resolution:-N_X_MIN) / M), Vector{Float64}(2 * pi * (N_X_MIN:resolution:N_X_MAX) / M))
# K_Y = vcat(Vector{Float64}(2 * pi * (-N_Y_MAX:resolution:-N_Y_MIN) / N), Vector{Float64}(2 * pi * (N_Y_MIN:resolution:N_Y_MAX) / N))


# pFTL = exp.(-im * [i * j for i in K_Y, j in 1:N])

# pFTR = exp.(-im * [i * j for i in 1:M, j in K_X])




# pftdat = abs.(pFTL * dat * pFTR)
# A = max(pftdat...)
# pftdat = pftdat*(1/A)
# downsampled = pftdat[1:10:size(pftdat)[1], 1:10:size(pftdat)[2]]
# downsampled = downsampled / max(downsampled...)
# Images.Gray.(downsampled)




# (ymax, xmax) = Tuple(findmax(pftdat)[2])
# (kx, ky) = (K_X[xmax], K_Y[ymax])


# ϕ = atan(ky, kx)
# ϕ = mod(ϕ + π / 2, π) - π / 2
# ϕ_deg = ϕ * (180 / pi)
# ϕ_vec[ind] = ϕ
# ϕ_vec_deg[ind] = ϕ_deg



