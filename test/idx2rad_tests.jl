# [test/idx2rad.jl]

using Test

@testset "Idx2Rad" begin
    import EarthAlbedo.idx2rad
    using MATLAB
    mat"""
        addpath('matlab_src');
    """   

    sy, sx = 288, 180

    N = 20
    is = range(0, sy, length = N)
    js = range(0, sx, length = N)
    for iIdx = 1:N
        for jIdx = 1:N

            i, j = Int(round(is[iIdx])), Int(round(js[jIdx]))

            mat"[$theta_mat, $eps_mat] = idx2rad($i, $j, $sy, $sx)"

            theta_jl, eps_jl = idx2rad(i, j, sy, sx)

            @test theta_mat == theta_jl 
            @test eps_mat == eps_jl

        end
    end
end

