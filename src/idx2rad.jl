# [scr/idx2rad.jl] 

""" 
Transforms from TOMS REFL Matrix indices to spherical coordinates (radians)

Arguments:
- i:  TOMS REFL Matrix latitude index of desired cell       |  Int
- j:  TOMS REFL Matrix longitude index of desired cell      |  Int 
- sy:  Number of latitude cells in grid                     |  Int 
- sx:  Number of longitude cells in grid                    |  Int 

Returns: (Tuple)
- (θ, ϵ):  Tuple of (Polar, Elevation) angle corresponding 
                to center of desired cell  (rad)             |  Tuple{Float, Float}
"""
function idx2rad(i::Int, j::Int, sy::Int, sx::Int)::Tuple{Float64, Float64}

    dx = 2 * pi / sx;
    dy = pi / sy;

    ϵ = pi - dy/2 - (i-1)*dy;
    θ = (j-1) * dx - pi + dx/2;

    return θ, ϵ
end