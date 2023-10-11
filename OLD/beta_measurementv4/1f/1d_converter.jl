include("/Users/alexjacoby/Documents/Research_Code/polarimeter/basic_functionality/1d_averager.jl")
import DelimitedFiles, Plots


n=23


small_fn_vec = vcat(reverse([string("-", i) for i in 1:20]),"0",[string("+", i) for i in 1:20])

small_fn = small_fn_vec[n]



cd("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv4/1f/1f_png/")


# vl = length(small_fn_vec)
# trough_1_vec = zeros(vl)
# trough_2_vec = zeros(vl)
# trough_3_vec = zeros(vl)

# deepest_trough_vec = zeros(vl)

# for (n,small_fn) in pairs(small_fn_vec)
#     raw_image = Images.Gray.(Images.load(string(small_fn, ".png")))[270:710, 700:1350]
#     dat = Float64.(raw_image)

#     correction = (Float64.(Images.Gray.(Images.load("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv4/background_low_intensity.png"))[270:710, 700:1350])) .^ (-1)
#     correction = (*(size(correction)...) / sum(correction)) * correction
#     dat = dat .* correction


#     index = DelimitedFiles.readdlm("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv4/1f/Estimated_Angles.txt")

#     angle = index[2,:][n]

#     (x,y) = convert_1d(dat, ϕ = angle)


#     min_vec = -log.(y)

#     deepest_trough_vec[n] = y[findmax(min_vec)[2]]

#     trough_1_vec[n] = y[100:700][findmax(min_vec[100:700])[2]]
#     trough_2_vec[n] = y[700:1300][findmax(min_vec[700:1300])[2]]
#     trough_3_vec[n] = y[1300:1900][findmax(min_vec[1300:1900])[2]]
# end

    
# Plots.plot(index[3,:], deepest_trough_vec)
# Plots.plot(index[3,:], trough_1_vec);
# Plots.plot!(index[3, :], trough_2_vec);
# Plots.plot!(index[3, :], trough_3_vec)

# Plots.plot(index[3,:], (trough_1_vec+trough_2_vec+trough_3_vec)/3)

# index[1,23]








raw_image = Images.Gray.(Images.load(string(small_fn, ".png")))[270:710, 700:1350];
dat = Float64.(raw_image)

correction = (Float64.(Images.Gray.(Images.load("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv4/background_low_intensity.png"))[270:710, 700:1350])) .^ (-1)
correction = (*(size(correction)...) / sum(correction)) * correction
dat = dat .* correction


index = DelimitedFiles.readdlm("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv4/1f/Estimated_Angles.txt")

angle = index[2, :][n]

(x, y) = convert_1d(dat, ϕ=angle)
x = x[100:1900]
y = y[100:1900]

Plots.plot(x,y)

model(r::Vector{Float64},parms::Vector{Float64}) = parms[1]*cos.(parms[3]*r .+parms[4]) .+ parms[2]


Plots.plot(x, model(x, [0.24, 0.24, 2.6, π - 0.3]))

parms0 = [0.24, 0.24, 2.6, π - 0.3]

x = Vector(0:2pi/(length(x)-1):2pi)
fit  = LsqFit.curve_fit(model,x, y, parms0)
fit.param
Plots.plot(x,y);
Plots.plot!(x, model(x,fit.param))
# Plots.plot!(x, [model(x_val,[0.24,0.24,2.6,π-0.3]) for x_val in x])