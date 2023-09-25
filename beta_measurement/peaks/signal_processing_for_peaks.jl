import DelimitedFiles, LinearAlgebra, Interpolations, FFTW, FourierAnalysis, Plots
function fermi(E::Number; β::Number, μ::Number)
    return 1 / (exp(β * (E - μ)) + 1)
end




number = 10



fn = string("beta_measurement/peaks/+", number, ".txt")



dat = DelimitedFiles.readdlm(fn)

x = dat[1, :] / max(abs.(dat[1, :])...)
y = dat[2, :]

Interpolations.linear_interpolation(x, y)
# Plots.plot(x,y)
analytic_y_fourier = FFTW.fft(FourierAnalysis.hilbert(y))

angle(analytic_y_fourier[6])
angle(analytic_y_fourier[11])




begin
    
    dat = DelimitedFiles.readdlm(fn)

    x = dat[1,:]/max(abs.(dat[1,:])...)
    y = dat[2,:]

    Interpolations.linear_interpolation(x,y)
    # Plots.plot(x,y)
    analytic_y_fourier = FFTW.fft(FourierAnalysis.hilbert(y))



    # Plots.plot!(abs.(FFTW.ifft(analytic_y_fourier .* second_filter))[300:1600], label = string(number))

    # second_filter = fermi.(1:length(y), β=-3.0, μ=9.0) .* fermi.(1:length(y), β=3.0, μ=15)




    fundamental_filter = fermi.(1:length(y), β=-3.0, μ=3.5) .* fermi.(1:length(y), β=3.0, μ=9)
    #plot the envelope to check against the 2f
    # Plots.plot(abs.(FFTW.ifft(analytic_y_fourier .* fundamental_filter)));
    #not yet clear to me whether this should be a sum normalization or a L^2 normalization
    first_envelope = LinearAlgebra.normalize((abs.(FFTW.ifft(analytic_y_fourier .* fundamental_filter))),1)*length(y)
    first_y = real.(FFTW.ifft(analytic_y_fourier .* fundamental_filter))






    second_filter = fermi.(1:length(y), β=-3.0, μ=9.) .* fermi.(1:length(y), β=3.0, μ=15)
    #plot the envelope to check against the 1f
    #Plots.plot(real.(FFTW.ifft(analytic_y_fourier .* second_filter)))

    second_envelope = LinearAlgebra.normalize(abs.(FFTW.ifft(analytic_y_fourier .* second_filter)),1)*length(y)
    second_y = real.(FFTW.ifft(analytic_y_fourier .* second_filter))
end
Plots.plot(x, first_y .* (first_envelope .^ (-1)) + second_y .* (first_envelope .^ (-1)))












