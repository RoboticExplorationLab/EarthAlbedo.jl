# [src/gridangle.jl]

""" 
Determines the angle between two cells given their indices 

Arguments: 
- i₁:  Latitude index of first cell in grid         |  Int
- j₁:  Longitude index of first cell in grid        |  Int 
- i₂:  Latitude index of second cell in grid        |  Int 
- j₂:  Longitude index of second cell in grid       |  Int
- sy:  Number of latitude cells in grid             |  Int 
- sx:  Number of longitude cells in grid            |  Int 

Outputs:
- ρ:   Angle between the two grid index pairs (rad) |  Float

"""
function gridangle(i₁, j₁, i₂, j₂, sy, sx)

    # Convert from index to angles
    θ₁, ϕ₁ = idx2rad(i₁, j₁, sy, sx);
    θ₂, ϕ₂ = idx2rad(i₂, j₂, sy, sx);

    # Determine the angular distance
    angle = sin(ϕ₁) * sin(ϕ₂) * cos(θ₁ - θ₂) + cos(ϕ₁) * cos(ϕ₂)
    
    # Correct numerical precision errors
    if (angle > 1.0) && (angle ≈ 1.0)
        angle = 1.0
    elseif (angle < -1.0) && (angle ≈ -1.0)
        angle = -1.0
    end

    return acos( angle )
end