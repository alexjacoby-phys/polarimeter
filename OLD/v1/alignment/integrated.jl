import Images, FileIO, LinearAlgebra
function rescale(x::Float64; α::Float64=1.0)
    return (exp(α * x) - 1) / (exp(α) - 1)
end


#fn_end = Base.prompt("Please Enter File Extension Under 'Data' Folder")


fn = "angle_benchmark/+4.png"
raw_image = Images.Gray.(Images.load(fn))[100:900, 300:1300]
dat = Float64.(raw_image)
(N, M) = size(dat)

#dat = [0.5*(1+cos((i+j)/100)) for i in 1:1080, j in 1:1440]


dat = dat - ones(Float64, size(dat)...) * sum(dat) / *(size(dat)...)

#This is the automatic version of this program. I don't recommend using it if your peaks are at all broad, or you can be very generous with δfringe if you do. It is advisable to do a coarse grained run with resolution at something like 0.01 and then a fine grained run at 0.001 or higher depending on the size of the peaks.



# guess_angle = -(pi / 180) * 45
# fringes = 5
# δfringe = 2
# resolution = 0.01

# N_X_MAX = abs((fringes+δfringe)*cos(guess_angle))
# N_Y_MAX = abs((fringes + δfringe) * sin(guess_angle))
# N_X_MIN = abs((fringes - δfringe) * cos(guess_angle))
# N_Y_MIN = abs((fringes - δfringe) * sin(guess_angle))


#1f
N_X_MAX = 4.1
N_X_MIN = 2.2
N_Y_MAX = 3.15
N_Y_MIN = 1.3
resolution = 0.001
# #2f
# N_X_MAX = 8.2
# N_X_MIN = 4.4
# N_Y_MAX = 6.3
# N_Y_MIN = 2.6
# resolution = 0.01


# K_X = Vector{Float64}(2 * pi * (-5:0.001:5) / M)
# K_Y = Vector{Float64}(2 * pi * (-5:0.001:5) / N)


K_X = vcat(Vector{Float64}(2 * pi * (-N_X_MAX:0.001:-N_X_MIN) / M), Vector{Float64}(2 * pi * (N_X_MIN:0.001:N_X_MAX) / M))
K_Y = vcat(Vector{Float64}(2 * pi * (-N_Y_MAX:0.001:-N_Y_MIN) / N), Vector{Float64}(2 * pi * (N_Y_MIN:0.001:N_Y_MAX) / N))


pFTL = exp.(-im * [i * j for i in K_Y, j in 1:N])

pFTR = exp.(-im * [i * j for i in 1:M, j in K_X])




pftdat = abs.(pFTL * dat * pFTR)
A = max(pftdat...)
pftdat = pftdat*(1/A)
downsampled = pftdat[1:10:size(pftdat)[1], 1:10:size(pftdat)[2]]
downsampled = downsampled / max(downsampled...)
Images.Gray.(downsampled)




(ymax, xmax) = Tuple(findmax(pftdat)[2])
(kx, ky) = (K_X[xmax], K_Y[ymax])


ϕ = atan(ky, kx)
ϕ = mod(ϕ + π / 2, π) - π / 2
ϕ_deg = ϕ * (180 / pi)



import Plots


kx_peak = pftdat[ymax, :]
ky_peak = pftdat[:, xmax]


FWHMX_vec = abs.(pftdat[ymax, :] .- 0.5 * pftdat[ymax, xmax])
FWHMY_vec = abs.(pftdat[:, xmax] .- 0.5 * pftdat[ymax, xmax])



# Plots.plot(K_X, kx_peak)
# Plots.plot(K_Y,ky_peak)

# Plots.plot(K_X, FWHMX_vec)
# Plots.plot(K_Y, FWHMY_vec)




FWHMXL = findmin(FWHMX_vec[1:xmax])[2]
FWHMXG = findmin(FWHMX_vec[xmax:length(FWHMX_vec)])[2] + xmax - 1
#+xmax-1 for FWHM-G


FWHMYL = findmin(FWHMY_vec[1:ymax])[2]
FWHMYG = findmin(FWHMY_vec[ymax:length(FWHMY_vec)])[2] + ymax - 1


FWHMX = abs(K_X[FWHMXG] - K_X[FWHMXL])
FWHMY = abs(K_Y[FWHMYG] - K_Y[FWHMYL])






import ImageTransformations





ImageTransformations.imrotate(raw_image, -ϕ)








# this sets the number of paritions for the one dimensional sample. Must be integer valued
itp = Interpolations.interpolate(dat, BSpline(Cubic(Line(OnGrid()))))
# how many data points will be on the x axis of the final one dimensional graph
no_partitions = 1000
#how much of a pixel should separate different sample locations in the average along constant phase lines (this should be less than or equal to one. Probably less is better.)
pixel_fraction = 1.

n1 = [cos(ϕ), -sin(ϕ)] *(pixel_fraction)



function isin_im(P::Vector{Float64};xmin = 1, xmax = 1440, ymin = 1, ymax = 1080 )
    b1 = (P[1] ≥ ymin) && (P[1] ≤ ymax)
    b2 = (P[2] ≥ xmin) && (P[2] ≤ xmax)
    return (b1 && b2)
end


#block for partitioning the k parallel lines

function point_partitions(p1::Vector{Float64},p2::Vector{Float64}; partitions::Int64 = 1)

    partition_coefficients = 0:(1/partitions):1.0

    return  [weight * p1 + (1 - weight) * p2 for weight in partition_coefficients]
end

function point_partitions(p1::Vector{Number},p2::Vector{Number}; partitions::Int64 = 1)

    partition_coefficients = 0:(1/partitions):1.0

    return  [weight * p1 + (1 - weight) * p2 for weight in partition_coefficients]
end


# returns a tuple with the first element a vector of points starting at p1 and ending at p2 equally spaced so that the steps are of length ℓ/partitions where ℓ is the distance between the two points. The second element is the distance from p1 ordered according to the list of points
function pp_d(p1::Vector{Float64}, p2::Vector{Float64}; partitions::Int64=1)
    ℓ = LinearAlgebra.norm(p1 - p2)
    partition_coefficients = 0:(1/partitions):1.0
    return ([(1 - weight) * p1 + weight * p2 for weight in partition_coefficients], Vector{Float64
    }(0:(ℓ/partitions):ℓ))
end

function pp_d(p1::Vector{Number}, p2::Vector{Number}; partitions::Int64=1)
    ℓ = LinearAlgebra.norm(p1-p2)
    partition_coefficients = 0:(1/partitions):1.0
    return ([(1 - weight) * p1 + weight * p2 for weight in partition_coefficients],Vector{Float64
    }(0:(ℓ/partitions):ℓ))
end
#maybe check the distances here




function pts_in_avg_slice(pavg::Vector{Float64}; normal::Vector{Float64})
    average_slice = [pavg]
    p = pavg + normal
    while isin_im(p)
        append!(average_slice, p)
        p += normal
    end
    p = pavg - normal
    while isin_im(p)
        append!(average_slice, p)
        p -= normal
    end
    return average_slice
end


function take_avg_of_slice(pavg::Vector{Float64}; normal::Vector{Float64}, dat_itp, xmin=1, xmax=1440, ymin=1, ymax=1080)
    average_slice = Vector{Float64}([])
    p = pavg
    while isin_im(p,xmin=1, xmax=M, ymin=1, ymax=N)
        append!(average_slice, dat_itp[p...])
        p += normal
    end
    p = pavg - normal
    while isin_im(p, xmin=1, xmax=M, ymin=1, ymax=N)
        append!(average_slice, dat_itp[p...])
        p -= normal
    end
    return sum(average_slice)/length(average_slice)
end



function avgpoints(no_partitions,ϕ; LX::Float64 = 1440.,LY::Float64 = 1080.)
    default_orientation = ϕ > 0


    #once things seem to work please replace this with N and M in case we change camera size/resolution


    ℓx = 0.5*(LX+1)
    ℓy = 0.5*(LY+1)

    p0 = (ℓy,ℓx)

    ap10x = (cot(ϕ)*ℓx + ℓy -1 + tan(ϕ))/(cot(ϕ)+tan(ϕ))
    ap10y = tan(ϕ)*(ap10x-1)+1

    ap20x = (tan(ϕ)*LX -LY + cot(ϕ)* ℓx + ℓy )/(tan(ϕ)+cot(ϕ))
    ap20y = tan(ϕ)*(ap20x- LX) + LY

    bp10x = (tan(ϕ)+ cot(ϕ)* ℓx+ ℓy - LY)/(cot(ϕ)+tan(ϕ))
    bp10y = tan(ϕ)*(bp10x-1)+LY

    bp20x = (cot(ϕ)*ℓx + tan(ϕ)* LX + ℓy - 1)/(cot(ϕ)+tan(ϕ))
    bp20y = tan(ϕ)*(bp20x - LX) + 1


    ap10 = [ap10y, ap10x]
    ap20 = [ap20y, ap20x]

    bp10 = [bp10y, bp10x]
    bp20 = [bp20y, bp20x]

    c1 = [1.,1]
    c2 = [1.,LX]
    c3 = Vector{Float64}([LY,LX])
    c4 = [LY,1.]

    if default_orientation
        (avgpoints1,dists1) = pp_d(ap10,c1; partitions = no_partitions);
        (avgpoints2, minusdists2) = pp_d(ap20,c3; partitions = no_partitions);
        dists2 = - minusdists2;
    else
        (avgpoints1, dists1) = pp_d(bp10, c4; partitions = no_partitions);
        (avgpoints2, minusdists2) = pp_d(bp20, c2; partitions = no_partitions);
        dists2 = -minusdists2;
    end
    return ((avgpoints1,dists1),(avgpoints2,dists2))
end


LX = M
LY = N
((avgpoints1,dists1),(avgpoints2,dists2)) = avgpoints(no_partitions,ϕ; LX = Float64(M), LY = Float64(N))


plotdat1 = take_avg_of_slice.(avgpoints1; normal=n1, dat_itp=itp)
plotdat2 = take_avg_of_slice.(avgpoints2; normal=n1, dat_itp=itp)

dists2

x = vcat(reverse(dists2),dists1[2:length(dists1)])
y = vcat(reverse(plotdat2),plotdat1[2:length(plotdat1)])

Plots.plot(x,y,legend = :none)




import Plots

Plots.plot(dists1, plotdat1)
Plots.plot!(dists2,plotdat2)


Plots.savefig("")


import DelimitedFiles
DelimitedFiles.writedlm("+4.txt",[x,y])


#LinearAlgebra.dot((c1 - ap10),n1) #good to convince yourself these are orthogonal (for default orientation = \phi positive)
#LinearAlgebra.dot((c4 - bp10),n1) #good to convince yourself these are orthogonal (for non-default orientation = \phi negative)







#this will make a really good figure for the writeup
avgpoints1_x = [val[2] for val in avgpoints1]
avgpoints1_y = [val[1] for val in avgpoints1]
avgpoints2_x = [val[2] for val in avgpoints2]
avgpoints2_y = [val[1] for val in avgpoints2]

Plots.plot(Images.Gray.(dat));
Plots.plot!(avgpoints1_x, avgpoints1_y);
Plots.plot!(avgpoints2_x, avgpoints2_y, legend = :none)









