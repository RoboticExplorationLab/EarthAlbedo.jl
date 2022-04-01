# [test/cellarea_tests.jl]

using Test 

@testset "Cell Area" begin 
    import EarthAlbedo.cellarea  
    using MATLAB
    mat"""
        addpath('matlab_src')
    """

    @testset "Default Dimensions" begin 
        sy, sx = 288, 180 
        
        is = Int.(round.(range(1, sy, length = 50)))
        js = Int.(round.(range(1, sx, length = 50)))

        for k = 1:length(is)
            i, j = is[k], js[k]

            mat"$a_mat = cellarea($i, $j, $sy, $sx)" 
            a_jl = cellarea(i, j, sy, sx)

            @test a_mat ≈ a_jl
        end
    end

    @testset "Changing a few things" begin 
        sy, sx = 50, 300 
        
        is = Int.(round.(range(1, sy, length = 50)))
        js = Int.(round.(range(1, sx, length = 50)))

        for k = 1:length(is)
            i, j = is[k], js[k]

            mat"$a_mat = cellarea($i, $j, $sy, $sx)" 
            a_jl = cellarea(i, j, sy, sx)

            @test a_mat ≈ a_jl
        end
    end


end

