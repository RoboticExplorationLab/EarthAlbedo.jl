# [src/earthfov.jl]

""" 
  Determines which cells on the Earth's surface are in the field of view 
of an object in space (e.g., satellite, sun) using spherical coordinates.

Arguments:
- pos_sph: Vector from Earth to the object in question (e.g., satellite, sun) 
            in ECEF frame using spherical coordinates [θ, ϵ, r]                      |  [3,]
- sy:      Number of latitude cells Earth's surface is divided into                  |  Int 
- sx:      Number of longitude cells Earth's surface is divided into                 |  Int 
- Rₑ:      (Optional) Average radius of the Earth (m)                                |  Float            

Returns:
- fov:  Cells on Earth's surface that are in the field of view of the given object   |  BitMatrix [sy, sx]
"""
function earthfov(pos_sph::Vector{Float64}, sy::Int, sx::Int, Rₑ = 6371.01e3)::BitMatrix # Matrix{Int8}

    if pos_sph[3] < Rₑ  # LEO shortcut (?) -
        pos_sph[3] = pos_sph[3] + Rₑ
        @warn "Warning: radial distance is less than radius of the Earth. Adding in Earth's radius..."
    end

    # Small circle center 
    θ₀, ϕ₀ = pos_sph[1], pos_sph[2]
    sϕ₀, cϕ₀ = sincos(ϕ₀)  # Precompute to use in the for loop

    ρ = acos(Rₑ / pos_sph[3]) # FOV on Earth 

    fov = BitArray(undef, sy, sx) # zeros(Int8, sy, sx)
    for i = 1:sy
        for j = 1:sx
            θ, ϕ = idx2rad(i, j, sy, sx)
            rd = acos( sϕ₀ * sin(ϕ) * cos(θ₀ - θ) + cϕ₀ * cos(ϕ) ); # Radial Distance

            fov[i, j] = (rd ≤ ρ) ? 1 : 0      # 1 if in FOV, 0 otherwise
        end 
    end

    return fov
end

""" 
  Determines which cells on the Earth's surface are in the field of view 
of an object in space (e.g., satellite, sun) using spherical coordinates.
This version takes in and updates a field-of-view (FOV) matrix in place.

Arguments:
- fov:     Cells on Earth's surface (updated in place with cells that are in field   |  Bitmatrix [sy, sx]
            of view of provided object)
- pos_sph: Vector from Earth to the object in question (e.g., satellite, sun)        |  [3,]
            in ECEF frame using spherical coordinates [θ, ϵ, r]                      
- sy:      Number of latitude cells Earth's surface is divided into                  |  Int 
- sx:      Number of longitude cells Earth's surface is divided into                 |  Int 
- Rₑ:      (Optional) Average radius of the Earth (m)                                |  Float            

"""
function earthfov!(fov::BitMatrix, pos_sph::SVector{3, Float64}, sy::Int, sx::Int, Rₑ = 6371.01e3)::BitMatrix # Matrix{Int8}

    if pos_sph[3] < Rₑ  # LEO shortcut 
        pos_sph[3] = pos_sph[3] + Rₑ
        @info "Warning: radial distance is less than radius of the Earth. Adding in Earth's radius..."
    end

    # Small circle center 
    θ₀, ϕ₀ = pos_sph[1], pos_sph[2]
    sϕ₀, cϕ₀ = sincos(ϕ₀)  # Precompute to use in the for loop

    ρ = acos(Rₑ / pos_sph[3]) # FOV on Earth 

    # fov = BitArray(undef, sy, sx) # zeros(Int8, sy, sx)
    for i = 1:sy
        for j = 1:sx
            if fov[i, j]
                θ, ϕ = idx2rad(i, j, sy, sx)
                rd = acos( sϕ₀ * sin(ϕ) * cos(θ₀ - θ) + cϕ₀ * cos(ϕ) ); # Radial Distance

                fov[i, j] = (rd ≤ ρ) ? 1 : 0      # 1 if in FOV, 0 otherwise
            end
        end 
    end

    return fov
end
