using LinearAlgebra, DelimitedFiles, Images, ImageTransformations, FileIO, Interpolations, LsqFit




function isin_im(P::Vector{Float64}, xmin, xmax, ymin, ymax)
    b1 = (P[1] ≥ ymin) && (P[1] ≤ ymax)
    b2 = (P[2] ≥ xmin) && (P[2] ≤ xmax)
    return (b1 && b2)
end


#block for partitioning the k parallel lines

function point_partitions(p1::Vector{Float64}, p2::Vector{Float64}; partitions::Int64=1)

    partition_coefficients = 0:(1/partitions):1.0

    return [weight * p1 + (1 - weight) * p2 for weight in partition_coefficients]
end

function point_partitions(p1::Vector{Number}, p2::Vector{Number}; partitions::Int64=1)

    partition_coefficients = 0:(1/partitions):1.0

    return [weight * p1 + (1 - weight) * p2 for weight in partition_coefficients]
end


# returns a tuple with the first element a vector of points starting at p1 and ending at p2 equally spaced so that the steps are of length ℓ/partitions where ℓ is the distance between the two points. The second element is the distance from p1 ordered according to the list of points
function pp_d(p1::Vector{Float64}, p2::Vector{Float64}; partitions::Int64=1)
    ℓ = LinearAlgebra.norm(p1 - p2)
    partition_coefficients = 0:(1/partitions):1.0
    return ([(1 - weight) * p1 + weight * p2 for weight in partition_coefficients], Vector{Float64
    }(0:(ℓ/partitions):ℓ))
end

function pp_d(p1::Vector{Number}, p2::Vector{Number}; partitions::Int64=1)
    ℓ = LinearAlgebra.norm(p1 - p2)
    partition_coefficients = 0:(1/partitions):1.0
    return ([(1 - weight) * p1 + weight * p2 for weight in partition_coefficients], Vector{Float64
    }(0:(ℓ/partitions):ℓ))
end





function pts_in_avg_slice(pavg::Vector{Float64}; normal::Vector{Float64}, A, B, C, D)
    nplus = 0
    nminus = 0

    p = pavg + normal
    while isin_im(p, xmin=A, xmax=B, ymin=C, ymax=D)
        p += normal
        nplus += 1
    end

    p = pavg - normal
    while isin_im(p, xmin=A, xmax=B, ymin=C, ymax=D)
        p -= normal
        nminus -= 1
    end

    return [pavg + k * normal for k in nminus:nplus]
end



function slice_avg(avg_point::Vector{Float64}, n1::Vector{Float64}; interpolation_obj, N, M)

    avg_pts = pts_in_avg_slice(avg_point, normal=n1, A=1, B=M, C=1, D=N)

    slice_vals = [interpolation_obj[p...] for p in avg_pts]

    return sum(slice_vals) / length(slice_vals)
end




function take_avg_of_slice(pavg::Vector{Float64}; normal::Vector{Float64}, dat_itp, xmin, xmax, ymin, ymax)
    average_slice = Vector{Float64}([])
    p = pavg
    while isin_im(p, xmin, xmax, ymin, ymax)
        append!(average_slice, dat_itp[p...])
        p += normal
    end
    p = pavg - normal
    while isin_im(p, xmin, xmax, ymin, ymax)
        append!(average_slice, dat_itp[p...])
        p -= normal
    end
    return sum(average_slice) / length(average_slice)
end



function avgpoints(no_partitions, ϕ; LX::Float64=1440.0, LY::Float64=1080.0)
    default_orientation = ϕ > 0


    #once things seem to work please replace this with N and M in case we change camera size/resolution


    ℓx = 0.5 * (LX + 1)
    ℓy = 0.5 * (LY + 1)

    p0 = (ℓy, ℓx)

    ap10x = (cot(ϕ) * ℓx + ℓy - 1 + tan(ϕ)) / (cot(ϕ) + tan(ϕ))
    ap10y = tan(ϕ) * (ap10x - 1) + 1

    ap20x = (tan(ϕ) * LX - LY + cot(ϕ) * ℓx + ℓy) / (tan(ϕ) + cot(ϕ))
    ap20y = tan(ϕ) * (ap20x - LX) + LY

    bp10x = (tan(ϕ) + cot(ϕ) * ℓx + ℓy - LY) / (cot(ϕ) + tan(ϕ))
    bp10y = tan(ϕ) * (bp10x - 1) + LY

    bp20x = (cot(ϕ) * ℓx + tan(ϕ) * LX + ℓy - 1) / (cot(ϕ) + tan(ϕ))
    bp20y = tan(ϕ) * (bp20x - LX) + 1


    ap10 = [ap10y, ap10x]
    ap20 = [ap20y, ap20x]

    bp10 = [bp10y, bp10x]
    bp20 = [bp20y, bp20x]

    c1 = [1.0, 1]
    c2 = [1.0, LX]
    c3 = Vector{Float64}([LY, LX])
    c4 = [LY, 1.0]

    if default_orientation
        (avgpoints1, dists1) = pp_d(ap10, c1; partitions=no_partitions)
        (avgpoints2, minusdists2) = pp_d(ap20, c3; partitions=no_partitions)
        dists2 = -minusdists2
    else
        (avgpoints1, dists1) = pp_d(bp10, c4; partitions=no_partitions)
        (avgpoints2, minusdists2) = pp_d(bp20, c2; partitions=no_partitions)
        dists2 = -minusdists2
    end
    return ((avgpoints1, dists1), (avgpoints2, dists2))
end







# this is the main 1d converter function. It accepts and requires the following arguments: dat (picture in matrix form), ϕ (angle), LX, LY (collectively the number of x, and y pixels. Note that this appears backwards, i.e., y comes first), no_partitions (number of points along the propagation direction-- how many points in x axis divided by two), pixel_fraction (determines the number of points that are averaged along an equal phase line to get a single y axis point)
function convert_1d(dat::Matrix{Float64}; ϕ::Float64, no_partitions::Int64 = 1000, pixel_fraction= 1.0)
    ϕ = ϕ = mod(ϕ + π / 2, π) - π / 2
    n1 = [cos(ϕ), -sin(ϕ)] * (pixel_fraction)
    itp = Interpolations.interpolate(dat, BSpline(Cubic(Line(OnGrid()))))
    (N,M) = size(dat)

    ((avgpoints1, dists1), (avgpoints2, dists2)) = avgpoints(no_partitions, ϕ; LX=Float64(M), LY=Float64(N))


    plotdat1 = take_avg_of_slice.(avgpoints1; normal=n1, dat_itp=itp, xmin = 1, xmax = M, ymin =  1, ymax = N)

    plotdat2 = take_avg_of_slice.(avgpoints2; normal=n1, dat_itp=itp, xmin = 1, xmax = M, ymin = 1, ymax = N)



    x = vcat(reverse(dists2), dists1[2:length(dists1)])
    y = vcat(reverse(plotdat2), plotdat1[2:length(plotdat1)])

    return (x,y)
end


#DON'T USE THIS ONE. IT IS SLOWER. this is the main 1d converter function. It accepts and requires the following arguments: dat (picture in matrix form), ϕ (angle), LX, LY (collectively the number of x, and y pixels. Note that this appears backwards, i.e., y comes first), no_partitions (number of points along the propagation direction-- how many points in x axis divided by two), pixel_fraction (determines the number of points that are averaged along an equal phase line to get a single y axis point)
function convert_1dv2(dat::Matrix{Float64}; ϕ::Float64, no_partitions::Int64=1000, pixel_fraction=1.0)
    ϕ = ϕ = mod(ϕ + π / 2, π) - π / 2
    n1 = [cos(ϕ), -sin(ϕ)] * (pixel_fraction)
    itp = Interpolations.interpolate(dat, BSpline(Cubic(Line(OnGrid()))))
    (N, M) = size(dat)

    ((avgpoints1, dists1), (avgpoints2, dists2)) = avgpoints(no_partitions, ϕ; LX=Float64(M), LY=Float64(N))

    plotdat1 = [slice_avg(pt, n1, interpolation_obj=itp, N=N, M=M) for pt in avgpoints1]

    plotdat2 = [slice_avg(pt, n1, interpolation_obj=itp, N=N, M=M) for pt in avgpoints2]

    x = vcat(reverse(dists2), dists1[2:length(dists1)])
    y = vcat(reverse(plotdat2), plotdat1[2:length(plotdat1)])

    return (x, y)
end







