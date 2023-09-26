import DelimitedFiles, LinearAlgebra, FFTW, FourierAnalysis, Plots, LsqFit
using Interpolations
function fermi(E::Number; β::Number, μ::Number)
    return 1 / (exp(β * (E - μ)) + 1)
end





number = 7
fn = string("beta_measurement/troughs/+", number, ".txt")


dat = DelimitedFiles.readdlm(fn)




x = dat[1, :] / max(abs.(dat[1, :])...)
y = dat[2, :]




Interpolations.linear_interpolation(x, y)
# Plots.plot(x,y)
analytic_y_fourier = FFTW.fft(FourierAnalysis.hilbert(y))




fundamental_filter = fermi.(1:length(y), β=-10.0, μ=4.0) .* fermi.(1:length(y), β=10.0, μ=9)
#plot the envelope to check against the 2f
# Plots.plot(abs.(FFTW.ifft(analytic_y_fourier .* fundamental_filter)));
#not yet clear to me whether this should be a sum normalization or a L^2 normalization
first_envelope = (abs.(FFTW.ifft(analytic_y_fourier .* fundamental_filter)))
first_envelope = first_envelope *(length(y)/sum(abs.(first_envelope)))
first_y = FFTW.ifft(analytic_y_fourier .* fundamental_filter)
y1 = first_y .* (first_envelope .^ (-1))
#Plots.plot(abs.(analytic_y_fourier .* fundamental_filter)[1:15])




second_filter = fermi.(1:length(y), β=-10.0, μ=9.0) .* fermi.(1:length(y), β=10.0, μ=15)
#plot the envelope to check against the 1f
#Plots.plot(real.(FFTW.ifft(analytic_y_fourier .* second_filter)))

second_envelope = abs.(FFTW.ifft(analytic_y_fourier .* second_filter))
second_envelope = second_envelope * (length(y)/sum(abs.(second_envelope)))
second_y = FFTW.ifft(analytic_y_fourier .* second_filter)
y2 = second_y .* (second_envelope .^ (-1))
#Plots.plot(abs.(analytic_y_fourier .* second_filter)[1:15])

Plots.plot(x,real.(y1+y2));


#parms = [k,A,B,ψ,ϕ]
@. model(r, parms::Vector{Float64}) = parms[2] * cos(parms[1] * r + parms[4]) + parms[3] * cos(2 * parms[1]* r + parms[5])
parms0 = [5.5 * pi,0.2,0.1,0,0]



fitting = LsqFit.curve_fit(model,x,real.(y1+y2),parms0)
result = fitting.param


result[4]- result[5]/2


Plots.plot(x,model(x, fitting.param));
Plots.plot!(x, real.(y1+y2))













































Plots.plot(x, first_y .* (first_envelope .^ (-1)));
Plots.plot!(x, second_y .* (first_envelope .^ (-1)))











number = 6
fn = string("beta_measurement/troughs/+", number, ".txt")
begin
    
    dat = DelimitedFiles.readdlm(fn)

    x = dat[1,:]/max(abs.(dat[1,:])...)
    y = dat[2,:]

    Interpolations.linear_interpolation(x,y)
    # Plots.plot(x,y)
    analytic_y_fourier = FFTW.fft(FourierAnalysis.hilbert(y))



    # Plots.plot!(abs.(FFTW.ifft(analytic_y_fourier .* second_filter))[300:1600], label = string(number))

    # second_filter = fermi.(1:length(y), β=-3.0, μ=9.0) .* fermi.(1:length(y), β=3.0, μ=15)




    fundamental_filter = fermi.(1:length(y), β=-10.0, μ=4.) .* fermi.(1:length(y), β=10.0, μ=9)
    #plot the envelope to check against the 2f
    # Plots.plot(abs.(FFTW.ifft(analytic_y_fourier .* fundamental_filter)));
    #not yet clear to me whether this should be a sum normalization or a L^2 normalization
    first_envelope = LinearAlgebra.normalize((abs.(FFTW.ifft(analytic_y_fourier .* fundamental_filter))),1)*length(y)
    first_y = real.(FFTW.ifft(analytic_y_fourier .* fundamental_filter))






    second_filter = fermi.(1:length(y), β=-10.0, μ=9.) .* fermi.(1:length(y), β=10.0, μ=15)
    #plot the envelope to check against the 1f
    #Plots.plot(real.(FFTW.ifft(analytic_y_fourier .* second_filter)))

    second_envelope = LinearAlgebra.normalize(abs.(FFTW.ifft(analytic_y_fourier .* second_filter)),1)*length(y)
    second_y = real.(FFTW.ifft(analytic_y_fourier .* second_filter))
end
Plots.plot(x, first_y .* (first_envelope .^ (-1))  );
Plots.plot!(x, second_y .* (first_envelope .^ (-1)))

Plots.plot(first_y[850:1000]);
findmax(first_y[850:1000])


Plots.plot(second_y[850:1000]);
findmax(second_y[850:1000])






x = 0:0.0001:pi

env(z) = 1-0.5*(2z/pi - 1)^2

env.(x)

k = 3*2pi
phi = -pi/13
y = cos.(k*x) #+ 0.5cos.(2k*x .+ phi)


signal = env.(x) .* y

abs.(FourierAnalysis.hilbert(signal))

Plots.plot(x, abs.(FourierAnalysis.hilbert(signal)))
