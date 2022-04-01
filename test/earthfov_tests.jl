# [test/earthfov.jl]

using Test 

@testset "Earth FOV" begin 
    import EarthAlbedo.earthfov 
    using CoordinateTransformations
    using MATLAB 
    mat"""
        addpath('matlab_src')
    """

    @testset "Defaults" begin
        sy, sx = 288, 180 
        
        # Run through some preset
        positions = [6.5e6  6.5e6  6.5e6;   # Near the surface
                     1.0e9  1.0e9  1.0e9;   # Far away 
                     1.0e3  2.0e3 -1.0e3;   # Too close (inside Earth) -> Should warn and then add in Re
                     0      0      7e6;     # +Z
                     0      0     -8e6;     # -Z
                     0      9e6    0;       # +Y 
                     0     -7e6    0;       # -Y
                     7.5e6  0      0;       # +X
                    -8.2e6  0      0        # -X
                    ]

        for i = 1:size(positions, 1)
            sat_sph_t = SphericalFromCartesian()(positions[i, :])
            sat_sph   = [sat_sph_t.θ; (pi/2) - sat_sph_t.ϕ; sat_sph_t.r];  # Adjust order to match template code to (θ ϵ r), with ϵ = (π/2) - ϕ

            mat"$fov_mat = earthfov($sat_sph, [$sy, $sx])"  # MATLAB wants them as a tuple

            fov_jl = earthfov(sat_sph, sy, sx)
            
            @test fov_mat == fov_jl
        end

        # # Throw in some more tests
        N = 50
        positions = ones(N, 3) * 6.5e6  + 3e6 * rand(N, 3) # Make at least Earth's radius away 
        signs = randn(N, 3)
        signs[signs .> 0] .= 1 
        signs[signs .<= 0] .= -1
        positions = positions .* signs

        for i = 1:size(positions, 1)
            sat_sph_t = SphericalFromCartesian()(positions[i, :])
            sat_sph   = [sat_sph_t.θ; (pi/2) - sat_sph_t.ϕ; sat_sph_t.r];  # Adjust order to match template code to (θ ϵ r), with ϵ = (π/2) - ϕ

            mat"$fov_mat = earthfov($sat_sph, [$sy, $sx])"  # MATLAB wants them as a tuple

            fov_jl = earthfov(sat_sph, sy, sx)
            
            @test fov_mat == fov_jl
        end

    end

    @testset "Modified Units" begin
        sy, sx = 20, 32 

        N = 50
        positions = ones(N, 3) * 6.5e6  + 3e6 * rand(N, 3) # Make at least Earth's radius away 
        signs = randn(N, 3)
        signs[signs .> 0] .= 1 
        signs[signs .<= 0] .= -1
        positions = positions .* signs ./ 1000 # Convert to kilometers
        Rₑ = 6371.01   # make km

        for i = 1:size(positions, 1)
            sat_sph_t = SphericalFromCartesian()(positions[i, :])
            sat_sph   = [sat_sph_t.θ; (pi/2) - sat_sph_t.ϕ; sat_sph_t.r];  # Adjust order to match template code to (θ ϵ r), with ϵ = (π/2) - ϕ

            sat_sph_unconverted = [sat_sph[1], sat_sph[2], 1000 * sat_sph[3]] 
            mat"$fov_mat = earthfov($sat_sph_unconverted, [$sy, $sx])"  # MATLAB wants them as a tuple

            fov_jl = earthfov(sat_sph, sy, sx, Rₑ)
            
            @test fov_mat == fov_jl
        end

    end

end
