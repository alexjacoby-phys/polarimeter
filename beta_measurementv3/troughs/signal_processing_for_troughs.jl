import DelimitedFiles, LinearAlgebra, FFTW, FourierAnalysis, Plots, LsqFit
using Interpolations
function fermi(E::Number; β::Number, μ::Number)
    return 1 / (exp(β * (E - μ)) + 1)
end

small_fn_vec = vcat(reverse([string("-", i) for i in 1:20]), "0", [string("+", i) for i in 1:20])

samples = length(small_fn_vec)

parameter_bank1= zeros(samples,5)

parameter_bank2 = zeros(samples,5)

error_bank1 = zeros(samples,5)

error_bank2 = zeros(samples,5)


for (k,number) in  pairs(small_fn_vec)
    
    fn = string("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv3/troughs/", number, ".txt")


    dat = DelimitedFiles.readdlm(fn)




    x = dat[1, :] / max(abs.(dat[1, :])...)
    y = dat[2, :]

    x = x[500:1500]
    y = y[500:1500]


    analytic_y_fourier = FFTW.fft(FourierAnalysis.hilbert(y))




    fundamental_filter = fermi.(1:length(y), β=-10.0, μ=4.0) .* fermi.(1:length(y), β=10.0, μ=9)
    #plot the envelope to check against the 2f
    # Plots.plot(abs.(FFTW.ifft(analytic_y_fourier .* fundamental_filter)));
    #not yet clear to me whether this should be a sum normalization or a L^2 normalization
    first_envelope = (abs.(FFTW.ifft(analytic_y_fourier .* fundamental_filter)))
    first_envelope = first_envelope * (length(y) / sum(abs.(first_envelope)))
    first_y = FFTW.ifft(analytic_y_fourier .* fundamental_filter)
    y1 = first_y .* (first_envelope .^ (-1))
    #Plots.plot(abs.(analytic_y_fourier .* fundamental_filter)[1:15])




    second_filter = fermi.(1:length(y), β=-10.0, μ=9.0) .* fermi.(1:length(y), β=10.0, μ=15)
    #plot the envelope to check against the 1f
    #Plots.plot(real.(FFTW.ifft(analytic_y_fourier .* second_filter)))

    second_envelope = abs.(FFTW.ifft(analytic_y_fourier .* second_filter))
    second_envelope = second_envelope * (length(y) / sum(abs.(second_envelope)))
    second_y = FFTW.ifft(analytic_y_fourier .* second_filter)
    y2 = second_y .* (second_envelope .^ (-1))
    #Plots.plot(abs.(analytic_y_fourier .* second_filter)[1:15])

    Plots.plot(x,real.(y1+y2))



    #parms = [k,A,B,ψ,ϕ]
    @. model(r, parms::Vector{Float64}) = parms[2] * cos(parms[1] * r + parms[4]) + parms[3] * cos(2 * parms[1] * r + parms[5])
    parms0 = [5.5 * pi, 0.1, 0.1, 0.0, 0.0]



    fitting1 = LsqFit.curve_fit(model, x, real.(first_y + second_y), parms0)
    result1 = fitting1.param
    parameter_bank1[k,:] = result1



    fitting2 = LsqFit.curve_fit(model, x, real.(y1 + y2), parms0)
    result2 = fitting2.param
    parameter_bank2[k,:] = result2

    error_bank1[k,:] = LsqFit.estimate_errors(fitting1)
    error_bank2[k,:] = LsqFit.estimate_errors(fitting2)
end

cd("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv3/troughs")

ϕ_vec = Vector{Float64}(DelimitedFiles.readdlm("angles.txt")[:,1])


parameter_bank1[:,4]
parameter_bank1[:,5]
parameter_bank1[21,4]
parameter_bank1[21,5]

Plots.plot(ϕ_vec, parameter_bank1[:, 5] .- 2*parameter_bank1[:, 4])
Plots.plot!(ϕ_vec,error_bank1[:,4]);
Plots.plot!(ϕ_vec, error_bank1[:, 5])






Plots.plot(ϕ_vec, mod.(parameter_bank2[:, 4] .+ parameter_bank1[:, 5], 2pi))
Plots.plot!(ϕ_vec, error_bank1[:, 4])
Plots.plot!(ϕ_vec, error_bank1[:, 5])











ϕ_vec[6]








number = string("-",2)
fn = string("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv3/troughs/", number, ".txt")


dat = DelimitedFiles.readdlm(fn)




x = dat[1, :] / max(abs.(dat[1, :])...)
y = dat[2, :]
x = x[500:1500]
y = y[500:1500]

analytic_y_fourier = FFTW.fft(FourierAnalysis.hilbert(y))




fundamental_filter = fermi.(1:length(y), β=-10.0, μ=4.0) .* fermi.(1:length(y), β=10.0, μ=9)
#plot the envelope to check against the 2f
# Plots.plot(abs.(FFTW.ifft(analytic_y_fourier .* fundamental_filter)));
#not yet clear to me whether this should be a sum normalization or a L^2 normalization
first_envelope = (abs.(FFTW.ifft(analytic_y_fourier .* fundamental_filter)))
first_envelope = first_envelope * (length(y) / sum(abs.(first_envelope)))
first_y = FFTW.ifft(analytic_y_fourier .* fundamental_filter)
y1 = first_y .* (first_envelope .^ (-1))
#Plots.plot(abs.(analytic_y_fourier .* fundamental_filter)[1:15])




second_filter = fermi.(1:length(y), β=-10.0, μ=9.0) .* fermi.(1:length(y), β=10.0, μ=15)
#plot the envelope to check against the 1f
#Plots.plot(real.(FFTW.ifft(analytic_y_fourier .* second_filter)))

second_envelope = abs.(FFTW.ifft(analytic_y_fourier .* second_filter))
second_envelope = second_envelope * (length(y) / sum(abs.(second_envelope)))
second_y = FFTW.ifft(analytic_y_fourier .* second_filter)
y2 = second_y .* (second_envelope .^ (-1))
#Plots.plot(abs.(analytic_y_fourier .* second_filter)[1:15])

#Plots.plot(x,real.(y1+y2))



#parms = [k,A,B,ψ,ϕ]
@. model(r, parms::Vector{Float64}) = parms[2] * cos(parms[1] * r + parms[4]) + parms[3] * cos(2 * parms[1] * r + parms[5])
parms0 = [5.5 * pi, 0.1, 0.1, 0.0, 0.0]



fitting1 = LsqFit.curve_fit(model, x, real.(first_y + second_y), parms0)
result1 = fitting1.param
result1



fitting2 = LsqFit.curve_fit(model, x, real.(y1 + y2), parms0)
result2 = fitting2.param


LsqFit.estimate_errors(fitting1)
LsqFit.estimate_errors(fitting2)



Plots.plot(x, real.(first_y + second_y), linedwidth=2, label="Fourier Filtered", fontfamily=:Times);
Plots.plot!(x, model(x, fitting1.param), linewidth=2, label="Fitted")
#Plots.savefig("Fourier_Filtered.pdf")




Plots.plot(x, real.(y1 + y2), label="Filtered + Evelope", linewidth=2, fontfamily=:Times);
Plots.plot!(x, model(x, fitting2.param), label="Fitted", linewidth=2)

#Plots.savefig("FilterandEnv.pdf")


Plots.plot(x, first_envelope, linewidth=2, legend=:none, fontfamily=:Times)
Plots.plot!(title="Envelope Function", ylabel="(Length-Normalized)")
#Plots.savefig("Envelope.pdf")







