import LinearAlgebra, Images, FileIO
#note that this interpolator must be imported as "using" not the more conventional "import." Some annoying stack of function dependencies happening.
using Interpolations





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




@time itp = Interpolations.interpolate(dat, BSpline(Cubic(Line(OnGrid()))))
n1 = [cos(ϕ), -sin(ϕ)]


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




function take_avg_of_slice(pavg::Vector{Float64}; normal::Vector{Float64},dat_itp)
    average_slice = Vector{Float64}([])
    p = pavg
    while isin_im(p)
        append!(average_slice, dat_itp[p...])
        p += normal
    end
    p = pavg - normal
    while isin_im(p)
        append!(average_slice, dat_itp[p...])
        p -= normal
    end
    return sum(average_slice)/length(average_slice)
end



# this sets the number of paritions for the one dimensional sample. Must be integer valued
no_partitions = 1000



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






take_avg_of_slice(avgpoints1[1]; normal = n1,dat_itp = itp)


plotdat1 = take_avg_of_slice.(avgpoints1; normal=n1, dat_itp=itp)
plotdat2 = take_avg_of_slice.(avgpoints2; normal=n1, dat_itp=itp)


import Plots

Plots.plot(dists1, plotdat1)
Plots.plot!(dists2,plotdat2)

#LinearAlgebra.dot((c1 - ap10),n1) #good to convince yourself these are orthogonal (for default orientation = \phi positive)
#LinearAlgebra.dot((c4 - bp10),n1) #good to convince yourself these are orthogonal (for non-default orientation = \phi negative)












