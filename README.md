SatelliteToolboxPropagators.jl
==============================

[![CI](https://github.com/JuliaSpace/SatelliteToolboxPropagators.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/JuliaSpace/SatelliteToolboxPropagators.jl/actions/workflows/ci.yml)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

This package defines orbit propagators for the **SatelliteToolbox.jl** ecosystem.

## Supported propagators

We currently support the following propagators:

| **Symbol**     | **Description**                |
|:---------------|:-------------------------------|
| `Val(:J2)`     | J2 orbit propagator algorithm. |

## Installation

``` julia
julia> using Pkg
julia> Pkg.install("SatelliteToolboxPropagators")
```
