import Images, FileIO, LinearAlgebra
function rescale(x::Float64; α::Float64=1.0)
    return (exp(α * x) - 1) / (exp(α) - 1)
end


#fn_end = Base.prompt("Please Enter File Extension Under 'Data' Folder")



fn = string("data/", "2f_test.png")
raw_image = Images.Gray.(Images.load(fn));
dat = Float64.(raw_image)
(N, M) = size(dat)

#dat = [0.5*(1+cos((i+j)/100)) for i in 1:1080, j in 1:1440]


dat = dat - ones(Float64, size(dat)...) * sum(dat) / *(size(dat)...)



# K_X = vcat(-5.:0.001:-4.,4:0.001:5)*(2pi/M)
# K_Y = vcat(-5.0:0.001:-4., 4:0.001:5)*(2pi/N)



K_X = 2 * pi * (-5:0.001:5) / M
K_Y = 2 * pi * (-5:0.001:5) / N


pFTL = exp.(-im * [i * j for i in K_Y, j in 1:N])

pFTR = exp.(-im * [i * j for i in 1:M, j in K_X])




pftdat = abs.(pFTL * dat * pFTR)
# A = max(pftdat...)
# pftdat = pftdat*(1/A)

# Images.Gray.(pftdat)




(ymax, xmax) = Tuple(findmax(pftdat)[2])
(kx, ky) = (K_X[xmax], K_Y[ymax])


ϕ = atan(ky, kx)
ϕ = mod(ϕ+π/2,π) - π/2






# anglemat = [mod(atan(ky, kx), pi) for ky in K_Y, kx in K_X]
# pftdat = rescale.(pftdat / max(pftdat...), α=20.0)

# pftdat = pftdat/sum(pftdat)
# θ =   sum(pftdat .* anglemat)

# import ImageTransformations


# mod(θ - ϕ,π)* (180/π)



# ImageTransformations.imrotate(raw_image, -θ)
ImageTransformations.imrotate(raw_image, -ϕ)

# rotdat = Float64.(ImageTransformations.imrotate(raw_image, θ))




# kx * (1440 / (2 * pi))
# ky * (1080 / (2 * pi))




# Images.Gray.(pftdat / max(pftdat...))
Images.Gray.(angle.(pFTL * dat * pFTR))






