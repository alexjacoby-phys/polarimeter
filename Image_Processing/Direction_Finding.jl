import Images, FileIO, LinearAlgebra
function rescale(x::Float64; α::Float64=1.0)
    return (exp(α * x) - 1) / (exp(α) - 1)
end


#fn_end = Base.prompt("Please Enter File Extension Under 'Data' Folder")



fn = string("early_test_data/", "1f_test.png")
raw_image = Images.Gray.(Images.load(fn))
dat = Float64.(raw_image)
(N, M) = size(dat)

#dat = [0.5*(1+cos((i+j)/100)) for i in 1:1080, j in 1:1440]


dat = dat - ones(Float64, size(dat)...) * sum(dat) / *(size(dat)...)

#This is the automatic version of this program. I don't recommend using it if your peaks are at all broad, or you can be very generous with δfringe if you do. It is advisable to do a coarse grained run with resolution at something like 0.01 and then a fine grained run at 0.001 or higher depending on the size of the peaks.



guess_angle = (pi/180)* 105
fringes = 3
δfringe = 3
resolution = 0.001

N_X_MAX = abs((fringes+δfringe)*cos(guess_angle))
N_Y_MAX = abs((fringes + δfringe) * sin(guess_angle))
N_X_MIN = abs((fringes - δfringe) * cos(guess_angle))
N_Y_MIN = abs((fringes - δfringe) * sin(guess_angle))







# K_X = Vector{Float64}(2 * pi * (-5:0.001:5) / M)
# K_Y = Vector{Float64}(2 * pi * (-5:0.001:5) / N)


K_X = vcat(Vector{Float64}(2 * pi * (-N_X_MAX:0.001:-N_X_MIN) / M), Vector{Float64}(2 * pi * (N_X_MIN:0.001:N_X_MAX) / M))
K_Y = vcat(Vector{Float64}(2 * pi * (-N_Y_MAX:0.001:-N_Y_MIN) / N), Vector{Float64}(2 * pi * (N_Y_MIN:0.001:N_Y_MAX) / N))


pFTL = exp.(-im * [i * j for i in K_Y, j in 1:N])

pFTR = exp.(-im * [i * j for i in 1:M, j in K_X])




pftdat = abs.(pFTL * dat * pFTR)
# A = max(pftdat...)
# pftdat = pftdat*(1/A)
downsampled = pftdat[1:10:size(pftdat)[1],1:10:size(pftdat)[2]]
downsampled = downsampled/max(downsampled...)
Images.Gray.(downsampled)




(ymax, xmax) = Tuple(findmax(pftdat)[2])
(kx, ky) = (K_X[xmax], K_Y[ymax])


import Plots


kx_peak = pftdat[ymax, :] 
ky_peak = pftdat[:,xmax] 


FWHMX_vec = abs.(pftdat[ymax, :] .- 0.5 * pftdat[ymax, xmax])
FWHMY_vec = abs.(pftdat[:, xmax] .- 0.5 * pftdat[ymax, xmax])



# Plots.plot(K_X, kx_peak)
# Plots.plot(K_Y,ky_peak)

# Plots.plot(K_X, FWHMX_vec)
# Plots.plot(K_Y, FWHMY_vec)




FWHMXL = findmin(FWHMX_vec[1:xmax])[2]
FWHMXG = findmin(FWHMX_vec[xmax:length(FWHMX_vec)])[2]+xmax-1
#+xmax-1 for FWHM-G


FWHMYL = findmin(FWHMY_vec[1:ymax])[2]
FWHMYG = findmin(FWHMY_vec[ymax:length(FWHMY_vec)])[2] + ymax - 1


FWHMX = abs(K_X[FWHMXG] - K_X[FWHMXL])
FWHMY = abs(K_Y[FWHMYG]-K_Y[FWHMYL])



ϕ = atan(ky, kx)
ϕ_deg = ϕ * (180 / pi)
ϕ = mod(ϕ+π/2,π) - π/2


# import ImageTransformations

# ImageTransformations.imrotate(raw_image, -ϕ)





