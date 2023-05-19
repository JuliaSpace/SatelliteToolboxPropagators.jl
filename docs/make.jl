using Documenter
using SatelliteToolboxPropagators

makedocs(
    modules = [SatelliteToolboxPropagators],
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),
        canonical = "https://juliaspace.github.io/SatelliteToolboxPropagators.jl/stable/",
    ),
    sitename = "Satellite Toolbox Propagators",
    authors = "Ronan Arraes Jardim Chagas",
    pages = [
        "Home" => "index.md",
        "Quick Start" => "man/quick_start.md",
        "Usage" => "man/usage.md",
        "Library" => "lib/library.md",
    ],
)

deploydocs(
    repo = "github.com/JuliaSpace/SatelliteToolboxPropagators.jl.git",
    target = "build",
)
