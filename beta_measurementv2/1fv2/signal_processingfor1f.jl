import DelimitedFiles, LinearAlgebra, Interpolations, FFTW, FourierAnalysis, Plots
function fermi(E::Number; β::Number, μ::Number)
    return 1 / (exp(β * (E - μ)) + 1)
end




list1 = []
list2 = []
list3 = []
small_fn_vec = vcat(reverse([string("+", i) for i in 1:10]), ["0"], [string("-", i) for i in 1:30])

for number in small_fn_vec
    fn = string("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv2/1fv2/", number, ".txt")
    begin

        dat = DelimitedFiles.readdlm(fn)

        x = dat[1, :] / max(abs.(dat[1, :])...)
        y = dat[2, :]

        Interpolations.linear_interpolation(x, y)
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




    #parms = [k,A,B,ψ,ϕ]
    @. model(r, parms::Vector{Float64}) = parms[2] * cos(parms[1] * r + parms[4]) + parms[3] * cos(2 * parms[1] * r + parms[5])
    parms0 = [5.5 * pi, 0.2, 0.01, 0.0, 0.0]



    fitting = LsqFit.curve_fit(model, x, real.(first_y + second_y), parms0)
    result1 = fitting.param


    result1[3]
    append!(list1, abs(result1[3]))

    fitting = LsqFit.curve_fit(model, x, real.(y1 + y2), parms0)
    result2 = fitting.param


    result2[3]
    append!(list2, abs(result2[3]))
end

list1
list2


Plots.plot( list1)
Plots.plot!( list2)





Plots.plot(x, model(x, fitting.param));
Plots.plot!(x, real.(first_y + second_y));
Plots.plot!(x, real.(y1 + y2))

Plots.plot(x, abs.(real.(first_y + second_y) - model(x, fitting.param)));
Plots.plot!(x, abs.(real.(y1 + y2) - model(x, fitting.param)))
















number = "-19"

fn = string("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv2/1fv2/", number, ".txt")
begin

    dat = DelimitedFiles.readdlm(fn)

    x = dat[1, :] / max(abs.(dat[1, :])...)
    y = dat[2, :]

    Interpolations.linear_interpolation(x, y)
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




#parms = [k,A,B,ψ,ϕ]
@. model(r, parms::Vector{Float64}) = parms[2] * cos(parms[1] * r + parms[4]) + parms[3] * cos(2 * parms[1] * r + parms[5])
parms0 = [5.5 * pi, 0.2, 0.01, 0.0, 0.0]



fitting = LsqFit.curve_fit(model, x, real.(first_y + second_y), parms0)
result1 = fitting.param


result1[3]
append!(list1, abs(result1[3]))

fitting = LsqFit.curve_fit(model, x, real.(y1 + y2), parms0)
result2 = fitting.param


result2[3]
append!(list2, abs(result2[3]))


Plots.plot(x, model(x, fitting.param));
Plots.plot!(x, real.(first_y + second_y));
Plots.plot!(x, real.(y1 + y2))

Plots.plot(x, abs.(real.(first_y + second_y) - model(x, fitting.param)));
Plots.plot!(x, abs.(real.(y1 + y2) - model(x, fitting.param)))