import Images, FileIO, LinearAlgebra
function rescale(x::Float64; α::Float64=1.0)
    return (exp(α * x) - 1) / (exp(α) - 1)
end

fn = string("data/", "2f_test.png")
raw_image = Images.Gray.(Images.load(fn));
dat = Float64.(raw_image)
(N, M) = size(dat)

#dat = [0.5*(1+cos((i+j)/100)) for i in 1:1080, j in 1:1440]


dat = dat - ones(Float64, size(dat)...) * sum(dat) / *(size(dat)...)



K_X = 2*pi*(-3:0.001:3)/M
K_Y = 2*pi*(-3:0.001:3)/N


pFTL = exp.(-im * [i * j for i in K_Y, j in 1:N])

pFTR = exp.(-im * [i * j for i in 1:M, j in K_X])




pftdat = abs.(pFTL * dat * pFTR)
pftdat = rescale.(pftdat/max(pftdat...),α=10.)


(ymax, xmax) = Tuple(findmax(pftdat)[2])
(kx,ky) = (K_X[xmax],K_Y[ymax])






θ = atan(ky,kx)


ImageTransformations.imrotate(raw_image,- θ)
rotdat = Float64.(ImageTransformations.imrotate(raw_image,θ))




kx * (1440/(2*pi))
ky * (1080/(2*pi))




Images.Gray.(pftdat/max(pftdat...))
Images.Gray.(angle.(pFTL * dat * pFTR))















L_X = 1:1440
L_Y = 1:1080
phi_rng = 
















pftdat = pftdat/sum(pftdat)





Images.Gray.(pftdat/max(pftdat...) )



anglemat = [ mod(atan(ky,kx),pi) for ky in K_Y, kx in K_X ]



sum(pftdat .* anglemat)*(180/pi)











[0.5*(cos(2*pi*(2i/1080-0.5(j/1440)))+1) for i in 1:1080, j in 1:1440]

x=1.2
Images.Gray.([0.5 * (cos(2 * pi *x* (2i / 1080 - 0.5(j / 1440))) + 1) for i in 1:1080, j in 1:1440])


Images.Gray.(([0.5 * (cos(2 * pi * x * (2i / 1080 - 0.5(j / 1440))) + 1) for i in 1:1080, j in 1:1440]+dat)*0.5)
