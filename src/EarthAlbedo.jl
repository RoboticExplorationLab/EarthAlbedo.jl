# [src/EarthAlbedo.jl]

module EarthAlbedo

    using LinearAlgebra 
    using SatelliteDynamics  # For GEOD/ECEF Conversions (may not work on Windows?)
    using CoordinateTransformations  # For spherical <-> Cartesian transformations
    using StaticArrays
    using JLD2

    include("rad2idx.jl")
    include("idx2rad.jl")
    include("gridangle.jl")
    include("earthfov.jl")
    include("cellarea.jl")

    export REFL, earth_albedo, load_refl

    """ REFL Structure for holding reflectivity data"""
    struct REFL
        data
        type
        start_time
        stop_time
    end

    """
      Determines the Earth's albedo at a satellite for a given set of reflectivity data. 
    First divides the Earth's surface into a set of cells (determined by the size of the refl data), 
    and then determines which cells are visible by both the satellite and the Sun. This is 
    then used to determine how much sunlight is reflected by the Earth towards the satellite.

    Arguments:
    - sat: Vector from center of Earth to satellite                                     |  SVector{3, Float64}
    - sun: Vector from center of Earth to Sun                                           |  SVector{3, Flaot64}
    - refl_data:  Matrix containing the averaged reflectivity data for each cell        |  Array{Float64, 2} 
    - Rₑ:   (Optional) Average radius of the Earth (default is 6371010 m)               |  Scalar 
    - AM₀:  (Optional) Solar irradiance  (default is 1366.9)                            |  Scalar

    Returns:
    - cell_albedos: Matrix with each element corresponding to the albedo coming from 
                        the corresponding cell on the surface of the earth            | [num_lat x num_lon]

    """
    function earth_albedo(sat::SVector{3, Float64}, sun::SVector{3, Float64}, refl_data::Array{Float64, 2}; Rₑ = 6371.01e3, AM₀ = 1366.9)   

        num_lat, num_lon = size(refl_data);                                            
        
        # Verify no satellite collision has occured 
        if norm(sat) < Rₑ                                                              
            error("albedo.m: The satellite has crashed into Earth!");
        end
    
        # # Convert from Cartesian -> Spherical (r θ ϕ)
        sat_sph_t = SphericalFromCartesian()(sat);                                      
        sun_sph_t = SphericalFromCartesian()(sun);                      
    
        # Adjust order to match template code to (θ ϵ r), with ϵ = (π/2) - ϕ
        sat_sph = SVector{3, Float64}(sat_sph_t.θ, (pi/2) - sat_sph_t.ϕ, sat_sph_t.r);  
        sun_sph = SVector{3, Float64}(sun_sph_t.θ, (pi/2) - sun_sph_t.ϕ, sun_sph_t.r);
    
        # # REFL Indices 
        sat_i, sat_j = rad2idx(sat_sph[1], sat_sph[2], num_lat, num_lon);          
        sun_i, sun_j = rad2idx(sun_sph[1], sun_sph[2], num_lat, num_lon);
    
        # # SKIP GENERATING PLOTS FOR NOW
    
        cells = fill!(BitArray(undef, num_lat, num_lon), 1)
        earthfov!(cells, sat_sph, num_lat, num_lon)
        earthfov!(cells, sun_sph, num_lat, num_lon)
    
        cell_albedos = zeros(Float64, num_lat, num_lon)                                                                             
    
        for i = 1:num_lat
            for j = 1:num_lon 
                if cells[i, j]
                    cell_albedos[i, j] = albedo_per_cell(sat, refl_data[i, j], i, j, sun_i, sun_j, num_lat, num_lon; Rₑ = Rₑ, AM₀ = AM₀)
                end
            end
        end
        
        return cell_albedos 
    end
    
    """
      Internal function that evaluates the albedo reflected by every cell.    

    Arguments:
    - sat: Vector from center of Earth to satellite                                     |  SVector{3, Float64}
    - refl_data:  Matrix containing the averaged reflectivity data for current cell     |  Array{Float64, 2} 
    - cell_i:     Latitude index of cell 
    - cell_j:     Longitude index of cell
    - sun_i:      Latitude index of sun
    - sun_j:      Longitude index of sun
    - num_lat:    Number of latitude cells Earth's surface is divided into              |  Int 
    - num_lon:    Number of longitude cells Earth's surface is divided into             |  Int 
    - Rₑ:   (Optional) Average radius of the Earth (default is 6371010 m)               |  Scalar 
    - AM₀:  (Optional) Solar irradiance  (default is 1366.9)                            |  Scalar

    Returns:
    - P_out:  Cell albedo generated by current cell towards the satellite               |  Float64
    """
    function albedo_per_cell(sat::SVector{3, Float64}, refl_cell_data::Float64, cell_i::Int, cell_j::Int, sun_i::Int, sun_j::Int, num_lat::Int, num_lon::Int; Rₑ = 6371.01e3, AM₀ = 1366.9)::Float64
        
        ϕ_in = gridangle(cell_i, cell_j, sun_i, sun_j, num_lat, num_lon);                  
        ϕ_in = (ϕ_in > pi/2) ? pi/2 : ϕ_in 
    
        E_inc = AM₀ * cellarea(cell_i, cell_j, num_lat, num_lon) * cos(ϕ_in)               
    
        grid_theta, grid_phi = idx2rad(cell_i, cell_j, num_lat, num_lon)                        
    
        grid_spherical = Spherical(Rₑ, grid_theta, (pi/2) - grid_phi)                                             
        grid = CartesianFromSpherical()(grid_spherical);                                                    
    
        sat_dist = norm(sat - grid);  # Unit vector pointing from grid center to satellite      
    
        ϕ_out = acos( ((sat - grid)/sat_dist)' * grid / norm(grid) ); # Angle to sat from grid  
    
        P_out = E_inc * refl_cell_data * cos(ϕ_out) / (pi  * sat_dist^2);    
    
        return P_out
    end


    """
      Loads in reflectivity data and stores it in a REFL struct
    """
    function load_refl(path = "data/refl.jld2")
        temp = load(path)
    
        refl = REFL( temp["data"], temp["type"], temp["start_time"], temp["stop_time"])
    
        return refl
    end

end # module
