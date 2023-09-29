import DelimitedFiles, LinearAlgebra, Interpolations, FFTW, FourierAnalysis, Plots
function fermi(E::Number; β::Number, μ::Number)
    return 1 / (exp(β * (E - μ)) + 1)
end


Plots.plot();

for number in 4:30
    fn = string("//Volumes/AJACOBY/beta_measurement/1f/+",number,".txt")
    dat = DelimitedFiles.readdlm(fn)

    x = dat[1,:]/max(abs.(dat[1,:])...)
    y = dat[2,:]

    Interpolations.linear_interpolation(x,y)
    # Plots.plot(x,y)
    analytic_y_fourier = FFTW.fft(FourierAnalysis.hilbert(y))



    filter_vec = fermi.(1:length(y), β=-1, μ=3.5)

    high_filter = fermi.(1:length(y), β = 1, μ = 20)
    append!(DUMBX, number)
    append!(DUMBY,sum(abs.(FFTW.ifft(analytic_y_fourier .* second_filter))[300:1600]))
end



Plots.plot(DUMBX,DUMBY,seriestype = :scatter)



DUMBX = []
DUMBY = []









number = 30
fn = string("//Volumes/AJACOBY/beta_measurement/1f/+", number, ".txt")
dat = DelimitedFiles.readdlm(fn)

x = dat[1, :] / max(abs.(dat[1, :])...)
y = dat[2, :]

Interpolations.linear_interpolation(x, y)
# Plots.plot(x,y)
analytic_y_fourier = FFTW.fft(FourierAnalysis.hilbert(y))



filter_vec = fermi.(1:length(y), β=-1, μ=3.5)

high_filter = fermi.(1:length(y), β=1, μ=20)


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

Plots.plot(x, first_y .* (first_envelope .^ (-1)) + second_y .* (second_envelope .^ (-1)))












x = 0:0.0001:pi

env(z) = 1-0.5*(2z/pi - 1)^2

env.(x)

k = 3*2pi
phi = -pi/13
y = cos.(k*x) #+ 0.5cos.(2k*x .+ phi)


signal = env.(x) .* y

abs.(FourierAnalysis.hilbert(signal))

Plots.plot(x, abs.(FourierAnalysis.hilbert(signal)))
