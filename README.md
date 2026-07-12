<p align="center">
  <img src="./docs/src/assets/logo.png" width="150" title="SatelliteToolboxTransformations.jl"><br>
  <small><i>This package is part of the <a href="https://github.com/JuliaSpace/SatelliteToolbox.jl">SatelliteToolbox.jl</a> ecosystem.</i></small>
</p>

# SatelliteToolboxPropagators.jl

[![CI](https://img.shields.io/github/actions/workflow/status/JuliaSpace/SatelliteToolboxPropagators.jl/ci.yml?style=flat-square&logo=githubactions&logoColor=white&labelColor=475569&label=CI)](https://github.com/JuliaSpace/SatelliteToolboxPropagators.jl/actions/workflows/ci.yml)
[![Codecov](https://img.shields.io/codecov/c/github/JuliaSpace/SatelliteToolboxPropagators.jl?token=WSVR7QYKOD&style=flat-square&logo=codecov&logoColor=white&labelColor=475569)](https://codecov.io/gh/JuliaSpace/SatelliteToolboxPropagators.jl)
[![docs-stable](https://img.shields.io/badge/docs-stable-16A34A?style=flat-square&logo=gitbook&logoColor=white&labelColor=475569)][docs-stable-url]
[![docs-dev](https://img.shields.io/badge/docs-dev-D97706?style=flat-square&logo=gitbook&logoColor=white&labelColor=475569)][docs-dev-url]
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495D1?style=flat-square&logo=julia&logoColor=white&labelColor=475569)](https://github.com/invenia/BlueStyle)
[![License](https://img.shields.io/github/license/JuliaSpace/SatelliteToolboxPropagators.jl?style=flat-square&logo=readme&logoColor=white&labelColor=475569&color=0284C7)](https://github.com/JuliaSpace/SatelliteToolboxPropagators.jl/blob/main/LICENSE.txt)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.11285127-DB2777?style=flat-square&logo=doi&logoColor=white&labelColor=475569)](https://zenodo.org/doi/10.5281/zenodo.11285127)

This packages contains orbit propagators for the **SatelliteToolbox.jl** ecosystem.

The current supported propagators are:

1. J2 analytical orbit propagator;
2. J2 osculating analytical orbit propagator;
3. J4 analytical orbit propagator;
4. J4 osculating analytical orbit propagator;
5. SGP4/SDP4 orbit propagator; and
6. Two body analytical orbit propagator.

## Installation

``` julia
julia> using Pkg
julia> Pkg.add("SatelliteToolboxPropagators")
```

## Documentation

For more information, see the [documentation][docs-stable-url].

[docs-dev-url]: https://juliaspace.github.io/SatelliteToolboxPropagators.jl/dev
[docs-stable-url]: https://juliaspace.github.io/SatelliteToolboxPropagators.jl/stable
