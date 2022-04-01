# [test/gridangle.jl]

using Test 

@testset "Gridangle" begin 
    import EarthAlbedo.gridangle 
    using MATLAB 
    mat"""
        addpath('matlab_src')
    """

    sy, sx = 288, 180 
    N = 500 

    i1s, j1s = rand(1:sy, N), rand(1:sx, N)
    i2s, j2s = rand(1:sy, N), rand(1:sx, N)
    
    for k = 1:N
        i1, j1 = i1s[k], j1s[k]
        i2, j2 = i2s[k], j2s[k]

        mat"$ρ_mat = gridangle($i1, $j1, $i2, $j2, $sy, $sx)"
        ρ_jl = gridangle(i1, j1, i2, j2, sy, sx)

        @test ρ_jl ≈ ρ_mat

    end
end


