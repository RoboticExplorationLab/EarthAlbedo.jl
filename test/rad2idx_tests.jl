# [test/rad2idx.jl]

using Test

@testset "Rad2Idx" begin
    import EarthAlbedo.rad2idx
    using JLD2

    @testset "Test 1" begin

        sy, sx = 288, 180

        N = 15
        θs = range(0, pi, length = N) # Cell center in spherical coordinates (radians)
        ϵs = range(0, 2 * pi, length = N)

        @load "test_files/rad2idx_test1.jld2" rad2idx_test1
        for i = 1:N
            for j = 1:N
                θ, ϵ = θs[i], ϵs[j]

                i_jl, j_jl = rad2idx(θ, ϵ, sy, sx)

                @test rad2idx_test1[i, j, 1] == i_jl 
                @test rad2idx_test1[i, j, 2] == j_jl      

            end
        end

    end;

    @testset "Test 2" begin

        sy, sx = 150, 3

        N = 15
        θs = range(0, pi, length = N) # Cell center in spherical coordinates (radians)
        ϵs = range(0, 2 * pi, length = N)

        @load "test_files/rad2idx_Test2.jld2" rad2idx_test2
        for i = 1:N
            for j = 1:N
                θ, ϵ = θs[i], ϵs[j]

                i_jl, j_jl = rad2idx(θ, ϵ, sy, sx)

                @test rad2idx_test2[i, j, 1] == i_jl 
                @test rad2idx_test2[i, j, 2] == j_jl        

            end
        end
    end

end

