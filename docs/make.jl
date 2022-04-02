using Documenter 
using EarthAlbedo 

makedocs(
    sitename = "EarthAlbedo.jl",
    format   = Documenter.HTML(prettyurls = false),
    pages    = [
        "Introduction" => "index.md",
        "API" => "api.md"
    ]
)

deploydocs(
    repo = "github.com/RoboticExplorationLab/EarthAlbedo.jl.git",
    devbranch = "main"
)