# [test/earthfov.jl]

using Test 

@testset "Earth FOV" begin 
    import EarthAlbedo.earthfov 
    using CoordinateTransformations
    using Random
    using JLD2

    @testset "Defaults" begin
        sy, sx = 288, 180 
        
        # Run through some preset positions
        positions = [6.5e6  6.5e6  6.5e6;   # Near the surface
                     1.0e9  1.0e9  1.0e9;   # Far away 
                     1.0e3  2.0e3 -1.0e3;   # Too close (inside Earth) -> Should warn and then add in Re
                     0      0      8e7;     # +Z
                     0      0     -8e7;     # -Z
                     0      9e7    0;       # +Y 
                     0     -7e8    0;       # -Y
                     7.5e6  0      0;       # +X
                    -8.2e6  0      0        # -X
                    ]

        N = size(positions, 1)

        @load "test_files/earthfov_test.jld2" earthfov_test
        for i = 1:N
            sat_sph_t = SphericalFromCartesian()(positions[i, :])
            sat_sph   = [sat_sph_t.θ; (pi/2) - sat_sph_t.ϕ; sat_sph_t.r];  # Adjust order to match template code to (θ ϵ r), with ϵ = (π/2) - ϕ

            if N == 3  # Should throw a warning...
                fov_jl = @test_logs (:warn, "Warning: radial distance is less than radius of the Earth. Adding in Earth's radius...") earthfov(sat_sph, sy, sx)
                @test earthfov_test[i] == fov_jl

            else
                fov_jl = earthfov(sat_sph, sy, sx)
                @test earthfov_test[i] == fov_jl
            end
        end
    end 

    @testset "More tests" begin
        sy, sx = 288, 180 

        # # Throw in some more tests
        N = 50
        Random.seed!(5)
        positions = ones(N, 3) * 6.5e6  + 3e6 * rand(N, 3) # Make at least Earth's radius away 
        signs = randn(N, 3)
        signs[signs .> 0] .= 1 
        signs[signs .<= 0] .= -1
        positions = positions .* signs

        @load "test_files/earthfov_more_tests.jld2" earthfov_more_tests
        for i = 1:N
            sat_sph_t = SphericalFromCartesian()(positions[i, :])
            sat_sph   = [sat_sph_t.θ; (pi/2) - sat_sph_t.ϕ; sat_sph_t.r];  # Adjust order to match template code to (θ ϵ r), with ϵ = (π/2) - ϕ

            fov_jl = earthfov(sat_sph, sy, sx)
            
            @test earthfov_more_tests[i] == fov_jl
        end

    end

    @testset "Modified Units" begin
        sy, sx = 20, 32 

        N = 50
        Random.seed!(5)
        positions = ones(N, 3) * 6.5e6  + 3e6 * rand(N, 3) # Make at least Earth's radius away 
        signs = randn(N, 3)
        signs[signs .> 0] .= 1 
        signs[signs .<= 0] .= -1
        positions = positions .* signs ./ 1000 # Convert to kilometers
        Rₑ = 6371.01   # make km


        @load "test_files/earthfov_mod_units.jld2" earthfov_mod_units
        for i = 1:N
            sat_sph_t = SphericalFromCartesian()(positions[i, :])
            sat_sph   = [sat_sph_t.θ; (pi/2) - sat_sph_t.ϕ; sat_sph_t.r];  # Adjust order to match template code to (θ ϵ r), with ϵ = (π/2) - ϕ

            sat_sph_unconverted = [sat_sph[1], sat_sph[2], 1000 * sat_sph[3]] 

            fov_jl = earthfov(sat_sph, sy, sx; Rₑ = Rₑ) 
            
            @test earthfov_mod_units[i] == fov_jl
        end

    end

end;
