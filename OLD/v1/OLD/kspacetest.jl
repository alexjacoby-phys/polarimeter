import LinearAlgebra, Plots

L = 1:1000

k = 2*pi*2/length(L)
dat = cos.(k*L)


dat = dat - sum(dat)*ones(size(dat)...)/*(size(dat)...)


Plots.plot(dat)






krange = 2*pi*Vector{Float64}(-5:0.0001:5)/100
[i * j for i in krange, j in L]
ftmat = exp.(im*[i*j for i in krange, j in L])



Plots.plot(krange,abs.(ftmat*dat))

pftdat = abs.(ftmat * dat)
pftdat = rescale.(pftdat/max(pftdat...),α=100.)
pftdat = pftdat*(1/sum(pftdat))

Plots.plot(krange, pftdat)
LinearAlgebra.dot(pftdat,abs.(krange))





function rescale(x::Float64;α::Float64=1.)
    return (exp(α * x)-1)/(exp(α)-1)
end