# [test/earth_albedo_tests.jl]

using Test
using StaticArrays

@testset "Earth Albedo Tests" begin 
    import EarthAlbedo.earth_albedo
    import EarthAlbedo.load_refl 
    import EarthAlbedo.REFL
    using JLD2, StaticArrays

    refl = load_refl("../data/refl.jld2");
    refl_data = refl.data;

    @testset "Standard" begin 
        sat = SVector{3, Float64}(0, 0, 8e6)
        sun = SVector{3, Float64}(0, 1e9, 0)

        sat_mat = [0, 0, 8e6]
        sun_mat = [0, 1e9, 0]
        
        @load "test_files/earth_alb_std.jld2" cell_albs_mat
        cell_albs_jl = earth_albedo(sat, sun, refl_data);

        @test cell_albs_jl ≈ cell_albs_mat atol=1e-12
    end

    @testset "Additional Positions" begin
        # different sat, sun positions
        sats = [  0     0    1e7;
                  0     1e8  0;
                  5e7   0    0;
                 -1e6   1e7  1e8;
                 -1e7  -1e7 -1e7;
                  0     1e9  0;
                  0     1e9  0  ]

        suns = [  0     0    1e9;
                  0     0    1e9;
                  0     0    1e9;
                  1e10  1e5  0;
                  -2e8  3300 9e9;
                  0     1e9  0;
                  0    -1e9  0  ]

        N = size(sats, 1)
        
        @load "test_files/earth_alb_extra.jld2" cell_albs_mats

        for i = 1:N
            sat_mat, sun_mat = sats[i, :], suns[i, :]
            sat = SVector{3, Float64}(sat_mat)
            sun = SVector{3, Float64}(sun_mat)

            cell_albs_jl = earth_albedo(sat, sun, refl_data);

            @test cell_albs_jl ≈ cell_albs_mats[i] atol=1e-12
        end
    end
   

end
