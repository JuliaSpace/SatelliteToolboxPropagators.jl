SatelliteToolboxPropagators.jl
==============================

[![CI](https://github.com/JuliaSpace/SatelliteToolboxPropagators.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/JuliaSpace/SatelliteToolboxPropagators.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/JuliaSpace/SatelliteToolboxPropagators.jl/branch/main/graph/badge.svg?token=WSVR7QYKOD)](https://codecov.io/gh/JuliaSpace/SatelliteToolboxPropagators.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

This package defines orbit propagators for the **SatelliteToolbox.jl** ecosystem.

## Supported propagators

We currently support the following propagators:

| **Symbol**     | **Description**                  |
|:---------------|:---------------------------------|
| `Val(:J2)`     | J2 orbit propagator algorithm.   |
| `Val(:J4)`     | J4 orbit propagator algorithm.   |
| `Val(:SGP4)`   | SGP4 orbit propagator algorithm. |

## Installation

``` julia
julia> using Pkg
julia> Pkg.install("SatelliteToolboxPropagators")
```
