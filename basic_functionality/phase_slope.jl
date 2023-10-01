import Images, FileIO, LinearAlgebra, LsqFit, ImageTransformations, OffsetArrays



raw_image = Images.Gray.(Images.load("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv4/1f/-1.png"))[270:710, 700:1350];

dat = Float64.(raw_image)


correction = (Float64.(Images.Gray.(Images.load("/Users/alexjacoby/Documents/Research_Code/polarimeter/beta_measurementv4/FFiltered_Background_Low.png"))[270:710, 700:1350])) .^ (-1)
correction = (*(size(correction)...) / sum(correction)) * correction;
dat = dat .* correction
Images.Gray.(dat)
function rotation(angle::Float64)
    return [cos(angle) -sin(angle); sin(angle) cos(angle)]
end


θ = -0.166

(N, M) = size(raw_image)
LX = M-1
LY = N-1


(PX, PY) = ((LX*cos(θ)-LY*sin(θ))/2, (LY*cos(θ)+LX*sin(θ))/2)



AR = 0.6
ψ = acos(1 / sqrt(1 + AR^2))
xval = (PX*tan(-θ)+PY)/(tan(ψ)+abs(tan(θ)))
R = Int64(floor((PX*tan(-θ)+PY)/(tan(ψ)+abs(tan(θ)))))
S = Int64(floor(xval*AR))

rotated = OffsetArrays.centered(ImageTransformations.imrotate(raw_image, θ));
strip_length =25
partitions = (2R) ÷ strip_length


sample = Float64.(rotated[-S:S,-R:R])
y = dropdims(sum(sample, dims = (1))/size(sample)[1], dims = (1))
x_array = Vector{Float64}(0:1/(length(y)-1):1)
C = sum(y)/length(y)
A = 2*sum(abs.(y .- C))/length(y)
B = 0 
parms0 = [4*pi, A, B, C]
model(r, parms::Vector{Float64}) = parms[2] * cos.(parms[1] * r .+ parms[3]) .+ parms[4]


params = zeros(partitions, 4)

Plots.plot(x_array,y)


for k in 1:partitions 
    indx = ((k-1)*strip_length+1):(k*strip_length)
    fitting = LsqFit.curve_fit(model, x_array[indx], y[indx], parms0)
    params[k, :] = fitting.param
    yfit = [model(x, fitting.param) for x in x_array[indx]]
    Plots.plot!(x_array[indx], yfit, legend=:none, color=:black)
end
Plots.plot!()

Plots.plot(1:partitions,mod.(params[:,1],2*pi))