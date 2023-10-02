
import Meshes



function rotation(angle::Float64)
    return [cos(angle) -sin(angle); sin(angle) cos(angle)]
end


# three keyword arguments: dims = size(image), angle = rotation angle, AR = final x size/ final y size
function rotation_crop(;dims::Tuple{Int64,Int64}, angle::Float64, AR::Number)
    (N, M) = dims
    LX = M - 1
    LY = N - 1


    C1 = [LY / 2, LX / 2]
    C2 = [-LY / 2, LX / 2]
    C3 = [-LY / 2, -LX / 2]
    C4 = [LY / 2, -LX / 2]

    rotmat = rotation(-angle)
    C1r = Tuple(rotmat * C1)
    C2r = Tuple(rotmat * C2)
    C3r = Tuple(rotmat * C3)
    C4r = Tuple(rotmat * C4)

    E12 = Meshes.Segment(C1r, C2r)
    E23 = Meshes.Segment(C2r, C3r)
    # E34 = Meshes.Segment(C3r, C4r)
    # E41 = Meshes.Segment(C4r, C1r)

    LR = Meshes.Line(Point(0, 0), Point(1.0, AR))
    RL = Meshes.Line(Point(0, 0), Point(1.0, -AR))


    LR_Intersections = [get(Meshes.intersection(LR, E12)), get(Meshes.intersection(LR, E23))]
    RL_Intersections = [get(Meshes.intersection(RL, E12)), get(Meshes.intersection(RL, E23))]
    LR_Intersections = abs.(collect.(Meshes.coordinates.(Vector{Meshes.Point2}(LR_Intersections[LR_Intersections.!=nothing])))[1])
    RL_Intersections = abs.(collect.(Meshes.coordinates.(Vector{Meshes.Point2}(RL_Intersections[RL_Intersections.!=nothing])))[1])


    LR_norm = LinearAlgebra.norm(LR_Intersections)

    RL_norm = LinearAlgebra.norm(RL_Intersections)


    if LR_norm > RL_norm
        R = Int64(floor(RL_Intersections[1]))
        S = Int64(floor(RL_Intersections[2]))
    elseif LR_norm < RL_norm
        R = Int64(floor(LR_Intersections[1]))
        S = Int64(floor(LR_Intersections[2]))
    elseif LR_norm == RL_norm
        R = Int64(floor(LR_Intersections[1]))
        S = Int64(floor(LR_Intersections[2]))
    end
    return (R,S)
end




# (N, M) = size(raw_image)
# LX = M - 1
# LY = N - 1


# C1 = [LY / 2, LX / 2]
# C2 = [-LY / 2, LX / 2]
# C3 = [-LY / 2, -LX / 2]
# C4 = [LY / 2, -LX / 2]

# rotmat = rotation(-θ)
# C1r = Tuple(rotmat * C1)
# C2r = Tuple(rotmat * C2)
# C3r = Tuple(rotmat * C3)
# C4r = Tuple(rotmat * C4)

# E12 = Meshes.Segment(C1r, C2r)
# E23 = Meshes.Segment(C2r, C3r)
# E34 = Meshes.Segment(C3r, C4r)
# E41 = Meshes.Segment(C4r, C1r)

# LR = Meshes.Line(Point(0, 0), Point(1.0, AR))
# RL = Meshes.Line(Point(0, 0), Point(1.0, -AR))


# LR_Intersections = [get(Meshes.intersection(LR, E12)), get(Meshes.intersection(LR, E23)), get(Meshes.intersection(LR, E34)), get(Meshes.intersection(LR, E41))]
# RL_Intersections = [get(Meshes.intersection(RL, E12)), get(Meshes.intersection(RL, E23)), get(Meshes.intersection(RL, E34)), get(Meshes.intersection(RL, E41))]
# LR_Intersections = collect.(Meshes.coordinates.(Vector{Meshes.Point2}(LR_Intersections[LR_Intersections.!=nothing])))
# RL_Intersections = collect.(Meshes.coordinates.(Vector{Meshes.Point2}(RL_Intersections[RL_Intersections.!=nothing])))


# LinearAlgebra.norm.(LR_Intersections)
# LinearAlgebra.norm.(RL_Intersections)
