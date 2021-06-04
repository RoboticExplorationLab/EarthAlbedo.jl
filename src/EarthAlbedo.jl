module EarthAlbedo
# Calculation of albedo for given satellite, sun location, and reflectivity data
# 
# (Ported over from Matlab's Earth Albedo Toolbox)
# REFL data downloaded from https://bhanderi.dk/downloads/


using CoordinateTransformations # For Cartesian/Spherical conversions
using SatelliteDynamics         # For GEOD/ECEF Conversions
using LinearAlgebra


CONSTANTS = (EARTH_RADIUS = 6371.01e3,    # Average radius of the Earth, m       
             SUN_IRRADIANCE = 1366.9,     # Irradiance of the sun (AM_0), W/m^2 
)

export refl_struct
struct refl_struct
    """ REFL Structure for holding reflectivity data"""
    data
    type
    start_time
    stop_time
end

export albedo
function albedo(sat, sun, refl)
    """
        Determines the Earth's albedo at a satellite for a given set of reflectivity data.
        - Divides Earth into a set of cells
        - Determines which cells are visible by both the satellite and the Sun
        - Determines how much sunlight is reflected by the Earth towards the satellite 

        Arguements:
        - sat: vector from center of Earth to satellite, [m]                              |  [3,]
        - sun: vector from center of Earth to Sun, [m]                                    |  [3,]
        - refl: Refl struct containing the averaged reflectivity data for each cell       | Custom struct

        Returns:
        - cell_albedos: Matrix with each element corresponding to the albedo coming from 
                            the corresponding cell on the surface of the earth            | [num_lat x num_lon]
        - union: Matrix containing indicators to which cells are illuminated by the sun 
                    AND visible by the satellite                                          | [num_lat x num_lon]

    """

    # Adjust input dimenions to support both row and column vectors
    sat = sat[:];
    sun = sun[:];

    num_lat, num_lon = size(refl.data); 

    # Verify no satellite collision has occured 
    if norm(sat) < CONSTANTS.EARTH_RADIUS
        error("albedo.m: The satellite has crashed into Earth!");
    end

    # Convert from Cartesian -> Spherical (r θ ϕ)
    sat_sph_t = SphericalFromCartesian()(sat);
    sun_sph_t = SphericalFromCartesian()(sun);

    # Adjust order to match template code to (θ ϵ r), with ϵ = (π/2) - ϕ
    sat_sph = [sat_sph_t.θ; (pi/2) - sat_sph_t.ϕ; sat_sph_t.r];
    sun_sph = [sun_sph_t.θ; (pi/2) - sun_sph_t.ϕ; sun_sph_t.r];

    # REFL Indices 
    sat_i, sat_j = rad2idx(sat_sph[1], sat_sph[2], num_lat, num_lon);
    sun_i, sun_j = rad2idx(sun_sph[1], sun_sph[2], num_lat, num_lon);

    # SKIP GENERATING PLOTS FOR NOW

    satellite_fov_cells = earthfov(sat_sph, num_lat, num_lon)   # Cells in the satellite field-of-view
    sunlit_cells     =    earthfov(sun_sph, num_lat, num_lon)   # Cells lit by the sun
    union = (satellite_fov_cells .!= 0) .& (sunlit_cells .!= 0) # Cells lit by the sun AND visible by the satellite

    cell_albedos = zeros(num_lat, num_lon)
    for i = 1:num_lat 
        for j = 1:num_lon 
            if union[i,j] # If a cell is illuminated by the sun AND visible by the satellite...
                ϕ_in = gridangle(i, j, sun_i, sun_j, num_lat, num_lon); # angle of incident solar Irradiance

                if ϕ_in > pi/2 # Account for numerical inaccuracies
                    ϕ_in = pi/2
                end

                E_inc = CONSTANTS.SUN_IRRADIANCE * cellarea(i, j, num_lat, num_lon) * cos(ϕ_in) # Incident power  

                grid_theta, grid_phi = idx2rad(i, j, num_lat, num_lon); # Distance to satellite from grid center

                grid_spherical = Spherical(CONSTANTS.EARTH_RADIUS, grid_theta, (pi/2) - grid_phi)
                grid = CartesianFromSpherical()(grid_spherical);    

                sat_dist = norm(sat - grid);  # Unit vector pointing from grid center to satellite

                ϕ_out = acos( ((sat - grid)/sat_dist)' * grid / norm(grid) ); # Angle to sat from grid

                P_out = E_inc * refl.data[i, j] * cos(ϕ_out) / (pi * sat_dist^2);

                cell_albedos[i, j] = P_out;
            end
        end
    end

    return cell_albedos, union;
end

function rad2idx(θ, ϵ, sy, sx)
    """ 
        Transforms location (in radians) to TOMS REFL matrix indices.
        
        Arguments:
        - θ, ϵ: Location of cell center in spherical coordinates (radians)                          | Scalars 
        - sy, sx: Number of latitude, longitude cells that the surface is being divided into        | Scalars 

        Returns: 
        - i, j: TOMS REFL matrix indices of cell corresponding to given spherical coordinates       | Scalars 
    """
    dx = 2 * pi / sx; # 360* / number of longitude cells 
    dy = pi / sy;     # 180* / number of latitude cells 

    i = round( (pi - dy/2 - ϵ)/dy ) + 1;
    j = round( (θ + pi - dx/2 )/dx ) + 1;

    # Adjust so that 180/-180 is included in the interval 
    if i == 0
        i = 1;
    end
    if j == 0;
        j = 1;
    end

    return i, j
end

function idx2rad(i, j, sy, sx)
    """ 
        Transforms TOMS REFL matrix indices to radians.

        Arguments:
        - i, j: TOMS REFL matrix indices of desired cell                                         | Scalars
        - sy, sx: Number of latitude, longitude cells that the surface is being divided into     | Scalars

        Returns:
        - θ, ϵ: Location of cell center in spherical coordinates (radians)                       | Scalars
    """

    dx = 2 * pi / sx;
    dy = pi / sy;

    ϵ = pi - dy/2 - (i-1)*dy;
    θ = (j-1) * dx - pi + dx/2;

    return θ, ϵ
end

function earthfov(pos_sph, sy, sx)
    """ 
        Determines the field of view on earth using spherical coordinates.

        Arguments:
        - pos_sph: vector from Earth to the object in question (satellite, sun) 
                        in ECEF frame using spherical coordinates                                      | [3,]
        - sy, sx: Number of latitude and longitude cells, respectively                                 | Scalars
    
        Returns:
        - fov: Cells on Earth's surface that are in the field of view of the given object (sat or sun) | [sy x sx]
    """

    IN_FOV = 1;         # Indicator that cell is in field-of-view
    NOT_IN_FOV = 0;     # Indicator that cell is NOT in field-of-view

    if pos_sph[3] < CONSTANTS.EARTH_RADIUS  # LEO shortcut (?)
        pos_sph[3] = pos_sph[3] + CONSTANTS.EARTH_RADIUS
    end

    dx = 2 * pi / sx;  
    dy = pi / sy;

    # Small circle center 
    θ_0 = pos_sph[1];
    ϕ_0 = pos_sph[2];

    ρ = acos(CONSTANTS.EARTH_RADIUS / pos_sph[3]) # FOV on Earth 
    fov = zeros(sy, sx)
    for i = 1:sy
        for j = 1:sx
            θ, ϕ = idx2rad(i, j, sy, sx)
            rd = acos( sin(ϕ_0)*sin(ϕ)*cos(θ_0-θ) + cos(ϕ_0)*cos(ϕ) ); # Radial Distance
            
            if rd <= ρ 
                fov[i,j] = IN_FOV;
            else
                fov[i,j] = NOT_IN_FOV;
            end
        end 
    end

    return fov
end

function gridangle(i1, j1, i2, j2, sy, sx)
    """ 
        Calculate the angle between two grid index pairs 

        Arguments:
        - i1, j1: First grid index pair                                       | Scalars
        - i2, j2: Second grid index pair                                      | Scalars
        - sy, sx: Number of latitude and longitude cells, respectively        | Scalars

        Returns:
        - ρ: Angle between the two grid index pairs                           | Scalars
    """

    θ1, ϕ1 = idx2rad(i1, j1, sy, sx)
    θ2, ϕ2 = idx2rad(i2, j2, sy, sx);

    ρ = acos( sin(ϕ1)*sin(ϕ2)*cos(θ1 - θ2) + cos(ϕ1)*cos(ϕ2) );

    return ρ
end

function cellarea(i, j, sy, sx)
    """ 
        Calculate area of TOMS cell for use in albedo calculation.
            
        Arguments:
        - i, j: TOMS REFL matrix indices of desired cell                                      | Scalars
        - sy, sx: Number of latitude and longitude cells, respectively                        | Scalars
        
        Returns:
        - A: Area of cell                                                                     | Scalars
    """
    _d2r = pi / 180.0; # Standard degrees to radians conversion

    θ, ϕ = idx2rad(i, j, sy, sx) # Convert to angles (radians)
    dϕ = (180.0 / sy) * _d2r;
    dθ = (360.0 / sx) * _d2r;

    # Get diagonal points
    ϕ_max = ϕ + dϕ/2;
    ϕ_min = ϕ - dϕ/2;

    A = (CONSTANTS.EARTH_RADIUS^2) * dθ * (cos(ϕ_min) - cos(ϕ_max));

    return A
end


export get_albedo_cell_centers
function get_albedo_cell_centers(lat_step = 1, lon_step = 1.25)
    """
        Returns the cell centers for the grid covering the surface of the Earth in Cartesian ECEF, to be used in later estimations of Earth's albedo,
            by looping through each cell's LLA coordinate and converting to ECEF 

        Arguments:
        - lat_step: (Optional) The step size (in degrees) to take across the latitude domain. Defaults to 1*        | Scalar 
        - lon_step: (Optional) The step size (in degrees) to take across the longitude domain. Defaults to 1.25*    | Scalar

        Returns:
        - cells_ecef: Matrix containing [x,y,z] coordinate for each latitude, longitude point.
                        Of form [lat, lon, [x,y,z]]                                                                 | [num_lat x num_lon x 3]
    """
    alt = 0.0 # Assume all cells are on surface of earth
    num_lat = Int(round((180 - lat_step) / lat_step) + 1)
    num_lon = Int(round((360 - lon_step) / lon_step) + 1)

    lon_offset = lon_step + (360 - lon_step) / 2   # Centers at 0 (longitude: [1.25, 360] => [-179.375, 179.375])
    lat_offset = lat_step + (180 - lat_step) / 2   # Centers at 0 (latitude:  [1.00, 180] => [-89.5, 89.5])

    cells_ecef = zeros(num_lat, num_lon, 3) # Lat, Lon, [x,y,z]
    for lat = 1:num_lat 
        for lon = 1:num_lon
            geod = [(lon * lon_step - lon_offset), (lat * lat_step - lat_offset), alt]
            ecef = sGEODtoECEF(geod, use_degrees = true)

            cells_ecef[Int(lat), Int(lon), :] = ecef
        end
    end

    return cells_ecef 
end

export get_diode_albedo
function get_diode_albedo(albedo_matrix, surface_normal, sat_pos)
    """ 
        Estimates the effect of Earth's albedo on a specific photodiode (by using the surface normal of that diode)
            = cell_albedo * surface_normal^T * r_g, with r_g as a unit vector in the direction of the grid point on Earth

        Arguments:
        - albedo_matrix: Albedo values for each cell on the Earth's surface         | [num_lat x num_lon] 
        - surface_normal: Photodiode surface normal                                 | [3,]
        - sat_pos: Cartesian position of satellite                                  | [3,]

        Returns:
        - diode_albedo: Total effect of albedo on specified photodiode              | Scalar
    """    
    cell_albedos = zeros(size(albedo_matrix))

    rows, cols = size(albedo_matrix)
    diode_albedo = 0.0
    for r = 1:1:rows
        for c = 1:1:cols
            if albedo_matrix[r,c] != 0
                r_g = cell_centers_ecef[r,c,:] - sat_pos # Distance from satellite to cell center
                r_g = r_g / norm(r_g)  # Make unit

                cell_albedo = (albedo_matrix[r,c] * (surface_normal * r_g))[1]

                if cell_albedo > 0.0    # Can't be negative
                    diode_albedo = diode_albedo + cell_albedo 
                end
            end
        end
    end
    
    return diode_albedo
end


end
