import DelimitedFiles, LinearAlgebra, FFTW, FourierAnalysis, Plots, LsqFit
using Interpolations
function fermi(E::Number; β::Number, μ::Number)
    return 1 / (exp(β * (E - μ)) + 1)
end


cd("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv3-2/peaks")
small_fn_vec = vcat(reverse([string("-", i) for i in 1:30]), "0", [string("+", i) for i in 1:20])




data = DelimitedFiles.readdlm.([string(fn, ".txt") for fn in small_fn_vec])

n = 6


Plots.plot(data[n][1, :], data[n][2, :])


index = DelimitedFiles.readdlm("angle_estimate.txt")

index[:, 7]