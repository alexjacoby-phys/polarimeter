using LinearAlgebra, LsqFit

ϕ = 0.0

k= 2.4*2pi
f(x, y) = (cos(k * cos(ϕ) * x + k * sin(ϕ)y) + cos(2*k * cos(ϕ) * x + 2*k * sin(ϕ)y) + 1.0) * 0.5

x_array = Vector{Float64}(0:0.001:1)
y_array = Vector{Float64}(0:0.001:1)
dat = [f(x,y) for x in x_array, y in y_array]


Images.Gray.(dat)

y = dropdims(sum(dat, dims = (2))/size(dat)[2], dims = (2))
x = x_array

Plots.plot(x,y)

strip_length = 20
partitions = length(y) ÷ strip_length



C = sum(y) / length(y)
A = 2 * sum(abs.(y .- C)) / length(y)
B = 0
parms0 = [8 * pi, A, B, C]
model(r, parms::Vector{Float64}) = parms[2] * cos.(parms[1] * r .+ parms[3]) .+ parms[4]


fitting = LsqFit.curve_fit(model, x_array, y, parms0)
Plots.plot(x_array, [model(x,fitting.param) for x in x_array])
Plots.plot!(x_array,y, linestyle = :dash)



params = zeros(partitions, 4)

Plots.plot(x_array, y);

n = 1
for k in n:(partitions-n)
    indx = ((k-n)*strip_length+1):((k+n)*strip_length)
    unit = 1:2*n*strip_length
    fitting = LsqFit.curve_fit(model, x_array[unit], y[indx], parms0)
    params[k, :] = fitting.param
    yfit = [model(x, fitting.param) for x in x_array[indx]]
    Plots.plot!(x_array[indx], yfit, legend=:none, color=:black)
end

Plots.plot(1:partitions,params[:,3])


for (k, val) in pairs(params[:, 1])
    if val < 0
        params[k, 3] = -params[k, 3]
        params[k, 1] = -params[k, 1]
    end
end
for (k, val) in pairs(params[:, 2])
    if val < 0
        params[k, 3] = params[k, 3] + pi
        params[k, 2] = -params[k, 2]
    end
end





Plots.plot(n:(partitions-n), mod.(params[n:(partitions-n), 3] .+ π, 2π) .- π)
