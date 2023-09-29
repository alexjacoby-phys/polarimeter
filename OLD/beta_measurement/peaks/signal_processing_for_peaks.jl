import DelimitedFiles, LinearAlgebra, Interpolations, FFTW, FourierAnalysis, Plots
function fermi(E::Number; β::Number, μ::Number)
    return 1 / (exp(β * (E - μ)) + 1)
end





# number = 10
# fn = string("beta_measurement/peaks/-", number, ".txt")


# dat = DelimitedFiles.readdlm(fn)



# x = dat[1, :] / max(abs.(dat[1, :])...)
# y = dat[2, :]

# scale = 10


# y_scaled = zeros(length(y) * scale)
# y_scaled[1:length(y)] = y




# analytic_y_f = FFTW.fft(FourierAnalysis.hilbert(y_scaled))



# Plots.plot(abs.(analytic_y_f[1:150]));


# findmax(abs.(analytic_y_f))

# angle(analytic_y_f[51])
# angle.(analytic_y_f[101]) / 2  -pi







number = 8
fn = string("beta_measurement/peaks/-", number, ".txt")
begin
    
    dat = DelimitedFiles.readdlm(fn)

    x = dat[1,:]/max(abs.(dat[1,:])...)
    y = dat[2,:]

    Interpolations.linear_interpolation(x,y)
    # Plots.plot(x,y)
    analytic_y_fourier = FFTW.fft(FourierAnalysis.hilbert(y))



    # Plots.plot!(abs.(FFTW.ifft(analytic_y_fourier .* second_filter))[300:1600], label = string(number))

    # second_filter = fermi.(1:length(y), β=-3.0, μ=9.0) .* fermi.(1:length(y), β=3.0, μ=15)




    fundamental_filter = fermi.(1:length(y), β=-10.0, μ=4.0) .* fermi.(1:length(y), β=10.0, μ=9)
    #plot the envelope to check against the 2f
    # Plots.plot(abs.(FFTW.ifft(analytic_y_fourier .* fundamental_filter)));
    #not yet clear to me whether this should be a sum normalization or a L^2 normalization
    first_envelope = (abs.(FFTW.ifft(analytic_y_fourier .* fundamental_filter)))
    first_envelope = first_envelope * (length(y) / sum(abs.(first_envelope)))
    first_y = FFTW.ifft(analytic_y_fourier .* fundamental_filter)
    y1 = first_y .* (first_envelope .^ (-1))





    second_filter = fermi.(1:length(y), β=-10.0, μ=9.0) .* fermi.(1:length(y), β=10.0, μ=15)
    #plot the envelope to check against the 1f
    #Plots.plot(real.(FFTW.ifft(analytic_y_fourier .* second_filter)))

    second_envelope = abs.(FFTW.ifft(analytic_y_fourier .* second_filter))
    second_envelope = second_envelope * (length(y) / sum(abs.(second_envelope)))
    second_y = FFTW.ifft(analytic_y_fourier .* second_filter)
    y2 = second_y .* (second_envelope .^ (-1))
end


data = real.(first_y+second_y)
#parms = [k,A,B,ψ,ϕ]
@. model(r, parms::Vector{Float64}) = parms[2] * cos(parms[1] * r + parms[4]) + parms[3] * cos(2 * parms[1] * r + parms[5])
parms0 = [5.5 * pi, 0.2, 0.15, 0.5, pi+0.5]



fitting = LsqFit.curve_fit(model, x, data, parms0)
result = fitting.param


mod(abs(result[4] - result[5] / 2),pi)  - pi/2


Plots.plot(x, model(x, fitting.param));
# Plots.plot(x, model(x, parms0));
Plots.plot!(x, data)










