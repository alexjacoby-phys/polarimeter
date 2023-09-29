import Images, FileIO, LinearAlgebra


fn = string("data/", "1f_test.png")
raw_image = Images.Gray.(Images.load(fn));
dat = Float64.(raw_image)
dat = dat - ones(Float64, size(dat)...) * sum(dat) / *(size(dat)...)
(N, M) = size(dat)

FTL = exp.(-2 * π * im * [i * j / N for i in 1:N, j in 1:N]) * (1 / sqrt(N))

FTR = exp.(-2 * π * im * [i * j / M for i in 1:M, j in 1:M]) * (1 / sqrt(M))




Images.Gray.(LinearAlgebra.abs.(FTL*dat*FTR))


Images.ImageView(log.(LinearAlgebra.norm.(FTL * dat * FTR)))



ftimg = Images.Gray.(log.(LinearAlgebra.norm.(FTL*dat*FTR)))[1:100,1:100]








fn = string("data/","1f_test.png")
raw_image = Images.Gray.(Images.load(fn))
dat = Float64.(raw_image)[1:1079,1:1439]
(N,M) = size(dat)

K = (N-1)÷2 
P = (M-1)÷2

L_N = -K:K
L_M = -P:P

FTL = exp.(-2 * π * im * [i * j / N for i in L_N, j in L_N]) * (1 / sqrt(N))

FTR = exp.(-2 * π * im * [i * j / M for i in L_M, j in L_M]) * (1 / sqrt(M))


dat_ft = LinearAlgebra.norm.(FTL * dat * FTR)

Images.Gray.(abs.(dat_ft)[5:540,700:740])












FTL = exp.(-2 * π * im * [i * j / N for i in 1:N, j in 1:N]) * (1 / sqrt(N))

FTR = exp.(-2 * π * im * [i * j / M for i in 1:M, j in 1:M])*(1/sqrt(M))


dat_ft = LinearAlgebra.norm.(FTL * dat * FTR)

Images.Gray.(abs.(dat_ft))















Images.Gray.(dat_ft/max(dat_ft...))

LinearAlgebra.svdvals(dat)



import FFTW



a = -0.1
b = 0.1
simulated = [0.5*(cos(a*i+b*j)+1) for i in 1:N, j in 1:M]
Images.Gray.(simulated)

simulated_ft = abs.(FTL * simulated * FTR)







Images.Gray.(abs.(simulated_ft))




imgmat = [(abs(x)+ abs(y))^(-0.5) for x in 0:100, y in 0:100]



Images.Gray.(imgmat)

