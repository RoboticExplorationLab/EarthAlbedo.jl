# [test/earth_albedo_tests.jl]

using Test
using StaticArrays

@testset "Earth Albedo Tests" begin 
    import EarthAlbedo.earth_albedo
    import EarthAlbedo.load_refl 
    import EarthAlbedo.REFL
    using MATLAB, MAT
    using JLD2, StaticArrays
    mat"""
        addpath('matlab_src')
    """

    function load_mat_refl(path = "../data/refl.mat")
        temp = matread(path)
    
        refl = REFL( temp["data"], temp["type"], temp["start_time"], temp["stop_time"])
    
        return refl
    end
    refl = load_refl("../data/refl.jld2");
    refl_data = refl.data;

    refl_mat = load_mat_refl()

    @testset "Standard" begin 
        sat = SVector{3, Float64}(0, 0, 8e6)
        sun = SVector{3, Float64}(0, 1e9, 0)

        sat_mat = [0, 0, 8e6]
        sun_mat = [0, 1e9, 0]
        mat"$cell_albs_mat = albedo($sat_mat, $sun_mat, $refl_mat)"
        cell_albs_jl = earth_albedo(sat, sun, refl_data);

        @test cell_albs_jl ≈ cell_albs_mat atol=1e-12
    end

    @testset "Additional Positions" begin
        # different sat, sun positions
        # Adjust size of refls? 
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
        for i = 1:N
            sat_mat, sun_mat = sats[i, :], suns[i, :]
            sat = SVector{3, Float64}(sat_mat)
            sun = SVector{3, Float64}(sun_mat)

            mat"$cell_albs_mat = albedo($sat_mat, $sun_mat, $refl_mat)"
            cell_albs_jl = earth_albedo(sat, sun, refl_data);

            @test cell_albs_jl ≈ cell_albs_mat atol=1e-12

        end
    end
   

end
