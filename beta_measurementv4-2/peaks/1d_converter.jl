include("/Users/alexjacoby/Documents/Research_Code/polarimeter/basic_functionality/1d_averager.jl")
import DelimitedFiles, Plots, LsqFit



small_fn_vec = vcat(reverse([string("-", i) for i in 1:30]), "0", [string("+", i) for i in 1:20])





cd("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv4/peaks/peaks_png/")


# n=11
small_fn = small_fn_vec[n]

vl = length(small_fn_vec)

parms = zeros(vl,6)



for (n, small_fn) in pairs(small_fn_vec)
    raw_image = Images.Gray.(Images.load(string(small_fn, ".png")))[270:710, 700:1350];
    dat = Float64.(raw_image)

    correction = (Float64.(Images.Gray.(Images.load("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv4/background_low_intensity.png"))[270:710, 700:1350])) .^ (-1)
    correction = (*(size(correction)...) / sum(correction)) * correction
    dat = dat .* correction


    index = DelimitedFiles.readdlm("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv4/peaks/Estimated_Angles.txt")



    angle = index[2, :][n]

    (x, y) = convert_1d(dat, ϕ=angle)
    Plots.plot(x, y)





    x = π * x/max(x...)


    model(r::Vector{Float64}, parms::Vector{Float64} )  = parms[1]*cos.(parms[4]*r .+ parms[5]) .+ parms[2]*cos.(2*parms[4]*r .+ parms[6]) .+ parms[3]


    parms0 = [0.05, 0.05, 0.2, 3.0, 0, π]

    fit = LsqFit.curve_fit(model, x, y, parms0, lower = [0.,0.,0.,2.,-π,-π], upper = [0.5,0.5,0.5, 3.5,π,π])
    parms[n,:] = fit.param

    # Plots.plot(x,y);
    # Plots.plot!(x,model(x,fit.param))

end


index[3,:]
Plots.plot(index[3,:],parms[:,5])
Plots.plot!(index[3,:],parms[:,6])



Plots.plot( (parms[:, 6] - 2 * parms[:,5])[20:30],color = :black);


Plots.plot!((π*ones(vl))[20:30], color = :red)


index[1,22]


