# [src/cellarea.jl]

""" 
Calculate the area of a TOMS REFL Matrix cell for use in calculating Earth albedo.

Arguments:
- i:  TOMS REFL Matrix latitude index of desired cell  (0 < i ≤ sy)         | Int
- j:  TOMS REFL Matrix longitude index of desired cell  (0 < j ≤ sx)        | Int
- sy:  Number of latitude cells in grid                                     | Int 
- sx:  Number of longitude cells in grid                                    | Int 
- Rₑ:  (Optional) Average radius of the Earth, defaults to meters           | Scalar

Returns:
- A:  Area of desired cell                                                  | Scalar
"""
function cellarea(i::Int, j::Int, sy::Int, sx::Int, Rₑ = 6371.01e3)    # Average radius of the Earth, m   

    _d2r = pi / 180.0; # Standard degrees to radians conversion

    θ, ϕ = idx2rad(i, j, sy, sx) # Convert to angles (radians)
    dϕ = (180.0 / sy) * _d2r;
    dθ = (360.0 / sx) * _d2r;

    # Get diagonal points
    ϕ_max = ϕ + dϕ/2;
    ϕ_min = ϕ - dϕ/2;

    A = (Rₑ^2) * dθ * (cos(ϕ_min) - cos(ϕ_max));

    return A
end