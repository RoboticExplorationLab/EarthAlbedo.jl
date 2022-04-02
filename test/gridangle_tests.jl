# [test/gridangle.jl]

using Test 

@testset "Gridangle" begin 
    import EarthAlbedo.gridangle
    using Random 

    sy, sx = 288, 180 
    N = 150

    Random.seed!(10)
    i1s, j1s = rand(1:sy, N), rand(1:sx, N)
    i2s, j2s = rand(1:sy, N), rand(1:sx, N)
    
    @load "test_files/gridangle_test.jld2" gridangle_test
    for k = 1:N
        i1, j1 = i1s[k], j1s[k]
        i2, j2 = i2s[k], j2s[k]

        ρ_jl = gridangle(i1, j1, i2, j2, sy, sx)

        @test ρ_jl ≈ gridangle_test[k]

    end
end


