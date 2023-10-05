include("/Users/alexjacoby/Documents/Research_Code/polarimeter/basic_functionality/1d_averager.jl")
import DelimitedFiles, Plots, LsqFit



small_fn_vec = vcat(reverse([string("-", i) for i in 1:20]), "0", [string("+", i) for i in 1:20])





cd("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv4/troughs/troughs_png/")

#n=37
#small_fn = small_fn_vec[n]

vl = length(small_fn_vec)

parms = zeros(vl, 6)



for (n, small_fn) in pairs(small_fn_vec)
    raw_image = Images.Gray.(Images.load(string(small_fn, ".png")))[270:710, 700:1350];
    dat = Float64.(raw_image)




    correction = (Float64.(Images.Gray.(Images.load("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv4/background_low_intensity.png"))[270:710, 700:1350])) .^ (-1)
    correction = (*(size(correction)...) / sum(correction)) * correction
    dat = dat .* correction


    index = DelimitedFiles.readdlm("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv4/troughs/Estimated_Angles.txt")



    aangle = index[2, :][n]

    (x, y) = convert_1d(dat, ϕ=aangle)
    x = x[250:1750]
    y = y[250:1750]



    Plots.plot(x, y)





    x = π * x / max(x...) 


    model(r::Vector{Float64}, parms::Vector{Float64}) = parms[1] * cos.(parms[4] * r .+ parms[5]) .+ parms[2] * cos.(2 * parms[4] * r .+ parms[6]) .+ parms[3]


    parms0 = [0.04, 0.06, 0.31, 2.37,0,0]# 0.12899908028433657-π, 0.48438720137510083 ]





    fit = LsqFit.curve_fit(model, x, y, parms0, lower=[0.03, 0.04, 0.3, 2.35, -Inf, -Inf], upper=[ 0.06, 0.08, 0.35, 2.4, Inf, Inf])
    parms[n, :] = fit.param


    # Plots.plot(x, y);
    #Plots.plot!(x, model(x, parms0));
    # Plots.plot!(x, model(x, fit.param))
end



Plots.plot(index[3, :], parms[:, 5]);
Plots.plot!(index[3, :], parms[:, 6])







rel_p = mod.(parms[:, 6] - 2 * parms[:, 5] .- π, 2π) .+ π
Plots.plot(rel_p - 2π * ones(vl));
Plots.plot!(2π*ones(vl)[10:20])

findmin(abs.(rel_p - 2π * ones(vl)))

abs.(rel_p - 2π * ones(vl))[17]

index[3,:16]