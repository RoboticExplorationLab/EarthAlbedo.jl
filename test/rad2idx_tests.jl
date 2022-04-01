# [test/rad2idx.jl]

using Test

@testset "Rad2Idx" begin
    import EarthAlbedo.rad2idx
    using MATLAB
    mat"""
        addpath('matlab_src');
    """   

    @testset "Test 1" begin

        sy, sx = 288, 180

        N = 15
        θs = range(0, pi, length = N) # Cell center in spherical coordinates (radians)
        ϵs = range(0, 2 * pi, length = N)
        for i = 1:N
            for j = 1:N
                θ, ϵ = θs[i], ϵs[j]
                mat"[$i_mat, $j_mat] = rad2idx($θ, $ϵ, $sy, $sx)"

                i_jl, j_jl = rad2idx(θ, ϵ, sy, sx)

                @test i_mat == i_jl 
                @test j_mat == j_jl      

            end
        end

    end

    @testset "Test 2" begin

        sy, sx = 150, 3

        N = 15
        θs = range(0, pi, length = N) # Cell center in spherical coordinates (radians)
        ϵs = range(0, 2 * pi, length = N)
        for i = 1:N
            for j = 1:N
                θ, ϵ = θs[i], ϵs[j]
                mat"[$i_mat, $j_mat] = rad2idx($θ, $ϵ, $sy, $sx)"

                i_jl, j_jl = rad2idx(θ, ϵ, sy, sx)

                @test i_mat == i_jl 
                @test j_mat == j_jl      

            end
        end
    end

end

