# [test/idx2rad.jl]

using Test

@testset "Idx2Rad" begin
    import EarthAlbedo.idx2rad
    using JLD2

    sy, sx = 288, 180

    N = 20
    is = range(0, sy, length = N)
    js = range(0, sx, length = N)

    @load "test_files/idx2rad_test.jld2" idx2rad_test
    for iIdx = 1:N
        for jIdx = 1:N

            i, j = Int(round(is[iIdx])), Int(round(js[jIdx]))
        
            theta_jl, eps_jl = idx2rad(i, j, sy, sx)

            @test idx2rad_test[iIdx, jIdx, 1] == theta_jl 
            @test idx2rad_test[iIdx, jIdx, 2] == eps_jl

        end
    end
end

