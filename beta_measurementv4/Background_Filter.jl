using Images, FFTW, LinearAlgebra

function fermi(E::Number; β::Number, μ::Number)
    return 1 / (exp(β * (E - μ)) + 1)
end


raw_image = Images.Gray.(Images.load("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv4/background_low_intensity.png"));
dat = Float64.(raw_image)
FT_DAT = FFTW.fftshift(FFTW.fft(dat))


(N,M) = size(dat)
distances = [LinearAlgebra.norm([y-N/2,x-M/2]) for y in 1:N, x in 1:M]
distances = distances/max(distances...)

filter = fermi.(distances,β = 20., μ = 0.3)

Images.Gray.(filter)

Images.save("FFiltered_Background_Low.png",Images.Gray.(abs.(FFTW.ifft(FFTW.fftshift(FT_DAT .* filter)))))
Images.Gray.(abs.(FFTW.ifft(FFTW.fftshift(FT_DAT .* filter))))


