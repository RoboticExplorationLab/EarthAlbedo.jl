# [test/runtests.jl]

using EarthAlbedo 
using Test 

using MATLAB

# Test scripts 
include("rad2idx_tests.jl")
include("idx2rad_tests.jl")
include("gridangle_tests.jl")
include("earthfov_tests.jl")
include("cellarea_tests.jl")
include("earth_albedo_tests.jl")