import DelimitedFiles, LinearAlgebra, Images, FileIO


small_fn = 7

fn = string("beta_measurement/troughs/+", small_fn)

raw_image = Images.Gray.(Images.load(string(fn, ".png")))[100:900, 300:1300];
dat = Float64.(raw_image)
(N, M) = size(dat)

#dat = [0.5*(1+cos((i+j)/100)) for i in 1:1080, j in 1:1440]


dat = dat - ones(Float64, size(dat)...) * sum(dat) / *(size(dat)...)




#1f
N_X_MAX = 2.9
N_X_MIN = 0.75
N_Y_MAX = 4
N_Y_MIN = 2.1
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









import FourierAnalysis, FFTW
x = -pi:0.001:pi

k = 1/rand()
phi = rand()*2pi
f(z) =  cos(k*z) + 0.5* cos(2*k*z+phi)


Plots.plot(x,f.(x))
y = f.(x)
y2 = cos.(k*x)
y3 = cos.(2k*x .+ phi)

ya = FourierAnalysis.hilbert(y)
y2a = FourierAnalysis.hilbert(y2)
y3a = FourierAnalysis.hilbert(y3)


yf = FFTW.fft(FourierAnalysis.hilbert(y))
y2f = FFTW.fft(FourierAnalysis.hilbert(y2))
y3f = FFTW.fft(FourierAnalysis.hilbert(y3))


Plots.plot(angle.(ya))
Plots.plot(angle.(y2a));
Plots.plot!(angle.(y3a))



Plots.plot(angle.(y2f[1:10]))
Plots.plot(angle.(yf[1:15]))


Plots.plot(angle.(y_p_f[1:100]),seriestype = :scatter)






Plots.plot(abs.(y_p_f)[1:100])

analytic_y_f = FFTW.fft(FourierAnalysis.hilbert(f.(x)))


Plots.plot(abs.(analytic_y_f)[1:100])



angle(analytic_y_f[5])











import Images, FileIO


Images.Gray.(Images.load("/Users/alexjacoby/Documents/Research_Code/polarimeter/OLD/1faligned.tif")) * 64
Images.Gray.(Images.load("/Users/alexjacoby/Documents/Research_Code/polarimeter/OLD/backgroundnowedge.tif")) * 64