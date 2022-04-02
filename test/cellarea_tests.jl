# [test/cellarea_tests.jl]

using Test 

@testset "Cell Area" begin 
    import EarthAlbedo.cellarea  
    using JLD2

    @testset "Default Dimensions" begin 
        sy, sx = 288, 180 
        N = 50
        
        is = Int.(round.(range(1, sy, length = N)))
        js = Int.(round.(range(1, sx, length = N)))

        @load "test_files/cellarea_default_dim.jld2" cellarea_default_dim
        for k = 1:length(is)
            i, j = is[k], js[k]

            a_jl = cellarea(i, j, sy, sx)

            @test cellarea_default_dim[k] ≈ a_jl
        end
    end;

    @testset "Changing a few things" begin 
        sy, sx = 50, 300 
        N = 50
        
        is = Int.(round.(range(1, sy, length = N)))
        js = Int.(round.(range(1, sx, length = N)))

        @load "test_files/cellarea_changes.jld2" cellarea_changes
        for k = 1:length(is)
            i, j = is[k], js[k]
            a_jl = cellarea(i, j, sy, sx)

            @test cellarea_changes[k] ≈ a_jl
        end
    end;


end;

