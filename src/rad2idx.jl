# [scr/rad2idx.jl] 

""" 
Transforms location of a cell center in spherical coordinates (radians) to TOMS REFL Matrix indices.
    (NOTE that we are forcing a specific rounding style to match MATLAB's)

Arguments:
- θ:  Polar angle corresponding to center of desired cell (rad)         |  Scalar 
- ϵ:  Elevation angle corresponding to center of desired cell (rad)     |  Scalar
- sy:  Number of latitude cells in grid                                 |  Int 
- sx:  Number of longitude cells in grid                                |  Int 

Returns: 
- (i, j):  TOMS REFL Matrix (latitude, longitude) index of desired cell  |  Tuple{Int, Int}
"""
function rad2idx(θ, ϵ, sy::Int, sx::Int)::Tuple{Int, Int}

    dx = 2.0 * pi / sx;  # 360* / number of longitude cells 
    dy = pi / sy;        # 180* / number of latitude cells 

    # Always round up to match MATLAB (vs default Bankers rounding)
    i = Int(round( (pi - dy/2.0 - ϵ)/dy  , RoundNearestTiesAway) + 1);
    j = Int(round( (θ + pi - dx/2.0 )/dx , RoundNearestTiesAway) + 1);

    return i, j
end