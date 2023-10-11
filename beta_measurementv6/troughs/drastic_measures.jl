import Images, FileIO, LinearAlgebra, LsqFit, ImageTransformations, OffsetArrays, Meshes, Plots
include("/Users/alexjacoby/Documents/Research_Code/polarimeter/basic_functionality/Rotation_Rectangle.jl")
include("/Users/alexjacoby/Documents/Research_Code/polarimeter/basic_functionality/1d_averager.jl")



small_fn_vec = small_fn_vec = vcat(reverse([string("-", i) for i in 1:25]), "0", [string("+", i) for i in 1:25])
γ = 68.5 * (π / 180)

n = 26
ϕ = 0.
small_fn = small_fn_vec[n]













# begin
raw_imagepr = Images.Gray.(Images.load(string("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv6/troughs/troughs_png/", small_fn, ".png")))

raw_image = OffsetArrays.centered(ImageTransformations.imrotate(raw_imagepr, γ))[-50:380, -100:400]
dat = Float64.(raw_image)

correction_imagepr = Images.Gray.(Images.load("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv6/background/backgroundv3.png"))
correction_image = OffsetArrays.centered(ImageTransformations.imrotate(correction_imagepr, γ))[-50:380, -100:400]

correction = Float64.(correction_image) .^ (-1)
correction = (*(size(correction)...) / sum(correction)) * correction;
dat = dat .* correction
image = Images.Gray.(dat)







(N, M) = size(dat)

#dat = [0.5*(1+cos((i+j)/100)) for i in 1:1080, j in 1:1440]


dat = dat - ones(Float64, size(dat)...) * sum(dat) / *(size(dat)...)




N_X_MAX = 5.5
N_X_MIN = 3
N_Y_MAX = 1
N_Y_MIN = 0
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
γ = (γ - ϕ)

γ*(180/pi)





















import Images, FileIO, LinearAlgebra, LsqFit, ImageTransformations, OffsetArrays, Meshes, Plots
include("/Users/alexjacoby/Documents/Research_Code/polarimeter/basic_functionality/Rotation_Rectangle.jl")
include("/Users/alexjacoby/Documents/Research_Code/polarimeter/basic_functionality/1d_averager.jl")



small_fn_vec = small_fn_vec = vcat(reverse([string("-", i) for i in 1:25]), "0", [string("+", i) for i in 1:25])
γ = -66.5 * (π / 180)
ϕ = 0.0







n = 26

small_fn = small_fn_vec[n]






n = 26

small_fn = small_fn_vec[n]








# begin
raw_imagepr = Images.Gray.(Images.load(string("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv6/peaks/peaks_png/", small_fn, ".png")))

raw_image = OffsetArrays.centered(ImageTransformations.imrotate(raw_imagepr, γ))[-450:0, -250:250]
dat = Float64.(raw_image)

correction_imagepr = Images.Gray.(Images.load("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv6/background/backgroundv3.png"))
correction_image = OffsetArrays.centered(ImageTransformations.imrotate(correction_imagepr, γ))[-450:0, -250:250]

correction = Float64.(correction_image) .^ (-1)
correction = (*(size(correction)...) / sum(correction)) * correction;
dat = dat .* correction
image = Images.Gray.(dat)
# end









(N, M) = size(dat)

#dat = [0.5*(1+cos((i+j)/100)) for i in 1:1080, j in 1:1440]


dat = dat - ones(Float64, size(dat)...) * sum(dat) / *(size(dat)...)




N_X_MAX = 4.5
N_X_MIN = 3.6
N_Y_MAX = 0.1
N_Y_MIN = 0
resolution = 0.0001




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

γ = γ - ϕ

γ*(180/pi)



