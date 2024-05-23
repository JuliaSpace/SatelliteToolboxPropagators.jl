using Documenter
using SatelliteToolboxPropagators

makedocs(
    modules = [SatelliteToolboxPropagators],
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),
        canonical = "https://juliaspace.github.io/SatelliteToolboxPropagators.jl/stable/",
        size_threshold = 500 * 1024,
        size_threshold_warn = 250 * 1024,
    ),
    sitename = "SatelliteToolboxPropagators.jl",
    authors = "Ronan Arraes Jardim Chagas",
    pages = [
        "Home" => "index.md",
        "Quick Start" => "man/quick_start.md",
        "Usage" => "man/usage.md",
        "Propagators" => [
            "J2" => "man/propagators/j2.md",
            "J2 Osculating" => "man/propagators/j2osc.md",
            "J4" => "man/propagators/j4.md",
            "J4 Osculating" => "man/propagators/j4osc.md",
            "SGP4/SDP4" => "man/propagators/sgp4.md",
            "Two-Body" => "man/propagators/twobody.md"
        ],
        "Library" => "lib/library.md",
    ],
)

deploydocs(
    repo = "github.com/JuliaSpace/SatelliteToolboxPropagators.jl.git",
    target = "build",
)
