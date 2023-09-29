using LinearAlgebra, Plots
begin
    L = 1000
    L_RNG = 1:1000
    N_RNG = 0.:0.0001:10

    k_rng = Vector((2π / L) * (N_RNG))
    ft = [exp(im * k * x) for k in k_rng, x in L_RNG]
    ftc = conj.(ft)


    k_vec = Vector{Float64}([])
    phi_vec = Vector{Float64}([])
    δk_vec = Vector{Float64}([])
    δϕ_vec = Vector{Float64}([])
    δk_f_vec = Vector{Float64}([])
    δϕ_f_vec = Vector{Float64}([])
end


for i in 1:100
    datk = max(k_rng...)*rand(Float64)
    append!(k_vec,datk)



    datphi = (2π)*rand(Float64)
    #Comment on this: I think that while the original binning is 2^8, I think we should be able to beat it when we do the average. It should be very possible to get someting like float16 out of it.
    data = #=Float16.=#(cos.(datk *L_RNG .+ datphi))

    begin
        ft_dat = ft*data
        ftc_dat = ftc*data
        density = abs.(ft_dat)

        density = density/max(density...)
        peak_ind = findmax(density)[2][1]

        w_pre = abs.(density .- 0.99)
        w_upper_ind = findmin(w_pre[peak_ind:length(w_pre)])[2][1]+peak_ind -1
        w_lower_ind = findmin(w_pre[1:peak_ind])[2][1]


        k_ffitted = k_rng[peak_ind]
        δk_f = abs(datk - k_ffitted)
        append!(δk_f_vec,δk_f)
        
        phi_ffitted = mod(0.5 * (angle(ftc_dat[peak_ind]) - angle(ft_dat[peak_ind])), 2pi)
        δϕ_f = abs(phi_ffitted- datphi)
        append!(δϕ_f_vec,δϕ_f)
    end




    k_fit_rng = k_rng[w_lower_ind]:(10e-7):k_rng[w_upper_ind]
    phi_fit_rng = (phi_ffitted - 0.15):(10e-5):(phi_ffitted + 0.15)


    full_fit = [ sum((cos.(k * L_RNG .+ ϕ) - data) .^2 ) for ϕ in phi_fit_rng, k in k_fit_rng]

    fitindices = Tuple(findmin(full_fit)[2])

    k_fitted = k_fit_rng[fitindices[2]]
    phi_fitted = phi_fit_rng[fitindices[1]]

    δϕ =  abs(datphi - phi_fitted)
    δk = abs(datk - k_fitted)
    append!(δk_vec,δk)
    append!(δϕ_vec,δϕ)

    println(i)
end



(L/(2*π))






Plots.plot(k_vec, -log10.(δϕ_vec), seriestype=:scatter)

Plots.plot(k_vec,-log10.(δk_vec),seriestype = :scatter)



Plots.plot(k_vec, δϕ_f_vec, seriestype=:scatter, ylims=(0, 0.2))


Plots.plot(k_vec,-log10.(δk_f_vec),seriestype = :scatter)


Plots.plot(k_vec , fwhm_vec,seriestype = :scatter)



k_error_avg = Vector{Float64}([])
phi_error_avg = Vector{Float64}([])

for (ind,val) in pairs(k_vec)
    if val > 0.025
        append!(phi_error_avg,δϕ_vec[ind])
        append!(k_error_avg,δk_vec[ind])
    end
end

sum(k_error_avg)/length(k_error_avg)

sum(phi_error_avg)/length(phi_error_avg)



