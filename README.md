<p align="center">
  <img src="./docs/src/assets/logo.png" width="150" title="SatelliteToolboxTransformations.jl"><br>
  <small><i>This package is part of the <a href="https://github.com/JuliaSpace/SatelliteToolbox.jl">SatelliteToolbox.jl</a> ecosystem.</i></small>
</p>

SatelliteToolboxPropagators.jl
==============================

[![CI](https://github.com/JuliaSpace/SatelliteToolboxPropagators.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/JuliaSpace/SatelliteToolboxPropagators.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/JuliaSpace/SatelliteToolboxPropagators.jl/branch/main/graph/badge.svg?token=WSVR7QYKOD)](https://codecov.io/gh/JuliaSpace/SatelliteToolboxPropagators.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)][docs-stable-url]
[![](https://img.shields.io/badge/docs-dev-blue.svg)][docs-dev-url]
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

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
julia> Pkg.install("SatelliteToolboxPropagators")
```

## Documentation

For more information, see the [documentation][docs-stable-url].

[docs-dev-url]: https://juliaspace.github.io/SatelliteToolboxPropagators.jl/dev
[docs-stable-url]: https://juliaspace.github.io/SatelliteToolboxPropagators.jl/stable
