import LinearAlgebra, OffsetArrays, Images



raw_image = Images.Gray.(Images.load("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv4/1f/-1.png"))[280:710, 700:1350];
dat = Float64.(raw_image)
correction = (Float64.(Images.Gray.(Images.load("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv4/FFiltered_Background_Low.png"))[280:710, 700:1350])) .^(-1)
correction = (*(size(correction)...)/sum(correction)) * correction;
dat = dat .* correction

Images.Gray.(dat);

N_X_MIN = 1.5
N_X_MAX = 3.5
N_Y_MIN = 0
N_Y_MAX = 2
resolution = 0.001





raw_image = 100 .* Images.Gray.(Images.load("/Users/alexjacoby/Documents/Research_Code/polarimeter/1faligned.tif"))[180:800, 750:1200];

dat = Float64.(raw_image)


correction = (Float64.(Images.Gray.(Images.load("/Users/alexjacoby/Documents/Research_Code/polarimeter/backgroundnowedge.tif"))[180:800, 750:1200])) .^ (-1)
correction = (*(size(correction)...) / sum(correction)) * correction;
dat = dat .* correction

Images.Gray.(dat)
#lawrence 1faligned
N_X_MIN = 0
N_X_MAX = 1.5
N_Y_MIN = 1.5
N_Y_MAX = 3.5
resolution = 0.001


(N, M) = size(dat)
dat = dat - ones(Float64, size(dat)...) * sum(dat) / *(size(dat)...)

K_X = vcat(Vector{Float64}(2 * pi * (-N_X_MAX:resolution:-N_X_MIN) / M), Vector{Float64}(2 * pi * (N_X_MIN:resolution:N_X_MAX) / M))
K_Y = vcat(Vector{Float64}(2 * pi * (-N_Y_MAX:resolution:-N_Y_MIN) / N), Vector{Float64}(2 * pi * (N_Y_MIN:resolution:N_Y_MAX) / N))


pFTL = exp.(-im * [i * j for i in K_Y, j in 1:N])

pFTR = exp.(-im * [i * j for i in 1:M, j in K_X])




pftdat = abs.(pFTL * dat * pFTR)
A = max(pftdat...)
pftdat = pftdat * (1 / A)
downsample_factor = 10
downsampled = pftdat[1:downsample_factor:size(pftdat)[1], 1:downsample_factor:size(pftdat)[2]]
downsampled = downsampled / max(downsampled...)
Images.Gray.(downsampled)



(ymax, xmax) = Tuple(findmax(pftdat[:,:])[2])
(kx, ky) = (K_X[xmax], K_Y[ymax])


ϕ = atan(ky, kx)
ϕ = mod(ϕ + π / 2, π) - π / 2
ϕ_deg = ϕ * (180 / pi)
ImageTransformations.imrotate(raw_image,-ϕ)








function peak_finder(dat::Matrix{Float64}; N_X_MIN,N_X_MAX,N_Y_MIN,N_Y_MAX,resolution::Float64,downsample_factor::Int64 = 10)
    (N, M) = size(dat)
    dat = dat - ones(Float64, size(dat)...) * sum(dat) / *(size(dat)...)

    K_X = vcat(Vector{Float64}(2 * pi * (-N_X_MAX:resolution:-N_X_MIN) / M), Vector{Float64}(2 * pi * (N_X_MIN:resolution:N_X_MAX) / M))
    K_Y = vcat(Vector{Float64}(2 * pi * (-N_Y_MAX:resolution:-N_Y_MIN) / N), Vector{Float64}(2 * pi * (N_Y_MIN:resolution:N_Y_MAX) / N))


    pFTL = exp.(-im * [i * j for i in K_Y, j in 1:N])

    pFTR = exp.(-im * [i * j for i in 1:M, j in K_X])




    pftdat = abs.(pFTL * dat * pFTR)
    A = max(pftdat...)
    pftdat = pftdat * (1 / A)
    downsampled = pftdat[1:downsample_factor:size(pftdat)[1], 1:downsample_factor:size(pftdat)[2]]
    downsampled = downsampled / max(downsampled...)
    return Images.Gray.(downsampled)
end



η = ϕ


ϕ_res = 0.001
δϕ = 10*ϕ_res

η_vec = Vector{Float64}((ϕ-δϕ):ϕ_res:(ϕ+δϕ))
σ = zeros(length(η_vec))


(x, y) = convert_1d(dat; ϕ=(mod(η+pi,pi)-π/2))
Plots.plot(x,y);

for (k,η) in pairs(η_vec)
    (x,y) = convert_1d(dat; ϕ = mod(η+π, π) - π/2)
    y = y[500:1500]
    mean = sum(y)/length(y)
    σ_vec = (y .- mean) .^2
    σ[k] = sum(σ_vec)
end

#import Plots
Plots.plot(η_vec, σ)



ϕ = η_vec[findmin(σ)[2]]
ϕ_res = ϕ_res/10
δϕ = 15*ϕ_res
η_vec = Vector{Float64}((ϕ-δϕ):ϕ_res:(ϕ+δϕ))
σ = zeros(length(η_vec))




ImageTransformations.imrotate(raw_image,-ϕ)