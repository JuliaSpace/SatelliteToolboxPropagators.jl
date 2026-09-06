# SGP4/SDP4

```@meta
CurrentModule = SatelliteToolboxPropagators
```

```@setup sgp4
using SatelliteToolboxPropagators
```

This package implements the interface to the
[SGP4/SDP4](https://en.wikipedia.org/wiki/Simplified_perturbations_models) propagator
provided by
[**SatelliteToolboxSgp4.jl**](https://github.com/JuliaSpace/SatelliteToolboxSgp4.jl).

## Algorithm

The SGP4/SDP4 implementation was built using **[1, 2, 3]**.

## Initialization

We must initialize the SGP4/SDP4 propagator with a [two-line element set
(TLE)](https://en.wikipedia.org/wiki/Two-line_element_set) using the following function:

```julia
Propagators.init(Val(:SGP4), tle::TLE; kwargs...) -> OrbitPropagatorSgp4
```

which creates a SGP4/SDP4 propagator structure [`OrbitPropagatorSgp4`](@ref) with the `tle`.
The following keyword selects the gravitational constants for the propagation algorithm:

- `sgp4c::Sgp4Constants`: SGP4 orbit propagator constants (see `Sgp4Constants`).
    (**Default**: `SGP4C_WGS84`)

The package [**SatelliteToolboxSgp4.jl**] contains some pre-built constants for this
propagator:

| **SGP4/SDP4 Propagator Constants** | **Description**           | **Type**  |
|-----------------------------------:|:--------------------------|:----------|
|                      `SGP4C_WGS84` | Constants based on WGS-84 | `Float64` |
|                      `SGP4C_WGS72` | Constants based on WGS-72 | `Float64` |

!!! note

    The type used in the propagation will be the same as used to define the constants in the
    structure `sgp4c`. The constants in another number type can be obtained with the
    converting constructor, e.g. `Sgp4Constants{Float32}(SGP4C_WGS84)`.

!!! note

    The package
    [**SatelliteToolboxTle.jl**](https://github.com/JuliaSpace/SatelliteToolboxTle.jl)
    defines the type `TLE`, which is re-exported here. It contains some useful
    functionalities, such as TLE fetching from online services. For more information, refer
    to the [package
    documentation](https://juliaspace.github.io/SatelliteToolboxTle.jl/stable/).

```@repl sgp4
tle = tle"""
       AMAZONIA 1
       1 47699U 21015A   23083.68657856 -.00000044  10000-8  43000-4 0  9990
       2 47699  98.4304 162.1097 0001247 136.2017 223.9283 14.40814394108652"""

Propagators.init(Val(:SGP4), tle)
```

We also support initializing an SGP4/SDP4 propagator passing the TLE information in
individual terms:

```julia
Propagators.init(Val(:SGP4), epoch::Number, n₀::Number, e₀::Number, i₀::Number, Ω₀::Number, ω₀::Number, M₀::Number, bstar::Number; kwargs...) -> OrbitPropagatorSgp4
```

where:

- `epoch::Number`: Epoch of the orbital elements [Julian Day].
- `n₀::Number`: SGP type "mean" mean motion at epoch [rad/s].
- `e₀::Number`: "Mean" eccentricity at epoch.
- `i₀::Number`: "Mean" inclination at epoch [rad].
- `Ω₀::Number`: "Mean" longitude of the ascending node at epoch [rad].
- `ω₀::Number`: "Mean" argument of perigee at epoch [rad].
- `M₀::Number`: "Mean" mean anomaly at epoch [rad].
- `bstar::Number`: Drag parameter (B*).

The keywords `kwargs...` are the same as in the first version of this function.

Finally, the propagator can be initialized with an [Orbit Mean-Elements Message
(OMM)](https://public.ccsds.org/Pubs/502x0b3e1.pdf) whose mean element theory is SGP4:

```julia
Propagators.init(Val(:SGP4), omm::OrbitMeanElementsMessage; kwargs...) -> OrbitPropagatorSgp4
```

The mean motion is obtained from the field `MEAN_MOTION` or, if it is absent, from the
semi-major axis and the gravitational coefficient, and the drag term is set to 0 if the
field `BSTAR` is absent. The keywords `kwargs...` are the same as in the first version of
this function.

!!! note

    The package
    [**SatelliteToolboxOrbitDataMessages.jl**](https://github.com/JuliaSpace/SatelliteToolboxOrbitDataMessages.jl)
    defines the type `OrbitMeanElementsMessage` and the functions to read, parse, and fetch
    the messages, and it is re-exported here.

```@repl sgp4
omm = parse_omm(
    """
    <?xml version="1.0" encoding="utf-8"?>
    <ndm><omm id="CCSDS_OMM_VERS" version="3.0">
    <header><CREATION_DATE>2025-12-30T23:36:37</CREATION_DATE><ORIGINATOR>18 SPCS</ORIGINATOR></header>
    <body><segment>
    <metadata><OBJECT_NAME>AMAZONIA 1</OBJECT_NAME><OBJECT_ID>2021-015A</OBJECT_ID><CENTER_NAME>EARTH</CENTER_NAME><REF_FRAME>TEME</REF_FRAME><TIME_SYSTEM>UTC</TIME_SYSTEM><MEAN_ELEMENT_THEORY>SGP4</MEAN_ELEMENT_THEORY></metadata>
    <data>
    <meanElements><EPOCH>2025-12-30T18:12:04.533984</EPOCH><MEAN_MOTION>14.40772474</MEAN_MOTION><ECCENTRICITY>0.00011240</ECCENTRICITY><INCLINATION>98.3721</INCLINATION><RA_OF_ASC_NODE>75.0877</RA_OF_ASC_NODE><ARG_OF_PERICENTER>97.3772</ARG_OF_PERICENTER><MEAN_ANOMALY>262.7545</MEAN_ANOMALY></meanElements>
    <tleParameters><EPHEMERIS_TYPE>0</EPHEMERIS_TYPE><CLASSIFICATION_TYPE>U</CLASSIFICATION_TYPE><NORAD_CAT_ID>47699</NORAD_CAT_ID><ELEMENT_SET_NO>999</ELEMENT_SET_NO><REV_AT_EPOCH>25439</REV_AT_EPOCH><BSTAR>0.00015330000000</BSTAR><MEAN_MOTION_DOT>0.00000447</MEAN_MOTION_DOT><MEAN_MOTION_DDOT>0.0000000000000</MEAN_MOTION_DDOT></tleParameters>
    </data>
    </segment></body>
    </omm></ndm>
    """
)

Propagators.init(Val(:SGP4), omm)
```

## Fitting Mean Elements

We can use the function:

```julia
Propagators.fit_mean_elements([sink::Type, ]::Val{:SGP4}, vjd::AbstractVector{Tjd}, vr_teme::AbstractVector{Tv}, vv_teme::AbstractVector{Tv}; kwargs...) -> sink, SMatrix{7, 7, Float64}
```

to fit a set of SGP4 mean elements using the osculating elements represented by a set of
position vectors `vr_teme` [m] and a set of velocity vectors `vv_teme` [m / s] represented
in the True-Equator, Mean-Equinox reference frame (TEME) at instants in the array `vjd`
[Julian Day, UTC]. The mean elements are returned as an object of type `sink`, which can be
a `TLE` or an `OrbitMeanElementsMessage`. If `sink` is omitted, an `OrbitMeanElementsMessage`
is returned.

It returns the fitted mean elements and the final covariance matrix of the least-square
algorithm.

This algorithm was based on **[4]**.

!!! note

    This algorithm version will allocate a new SGP4 propagator with the default constants
    `SGP4C_WGS84`. If another set of constants are required, use the function
    [`Propagators.fit_mean_elements!`](@ref) instead.

The following keywords are available to configure the fitting process:

- `atol::Number`: Tolerance for the residual absolute value. If the residual is lower than
    `atol` at any iteration, the computation loop stops.
    (**Default**: 2e-4)
- `rtol::Number`: Tolerance for the relative difference between the residuals. If the
    relative difference between the residuals in two consecutive iterations is lower than
    `rtol`, the computation loop stops.
    (**Default**: 2e-4)
- `estimate_bstar::Bool`: If `true`, the algorithm will try to estimate the B* parameter.
    Otherwise, it will be set to 0 or to the value in the initial guess.
    (**Default**: `true`)
- `include_covariance::Bool`: If `true`, the covariance of the mean position and velocity
    obtained by the least-square algorithm is stored in the covariance matrix section of
    the output, represented in the TEME reference frame. It is only used when `sink` is
    `OrbitMeanElementsMessage`.
    (**Default**: `true`)
- `initial_guess::Union{Nothing, AbstractVector, TLE, OrbitMeanElementsMessage}`: Initial
    guess for the fitting process. If it is `nothing`, the algorithm will obtain an initial
    estimate from the osculating elements in `vr_teme` and `vv_teme`.
    (**Default**: `nothing`)
- `jacobian_method::AbstractJacobianMethod`: Method used to compute the Jacobian matrix. It
    can be `FiniteDiffJacobian()` for finite differences or `ForwardDiffJacobian()` for
    **ForwardDiff.jl** automatic differentiation.
    (**Default**: `FiniteDiffJacobian()`)
- `jacobian_perturbation::Number`: Initial state perturbation to compute the
    finite-difference when calculating the Jacobian matrix.
    (**Default**: 1e-3)
- `jacobian_perturbation_tol::Number`: Tolerance to accept the perturbation when calculating
    the Jacobian matrix. If the computed perturbation is lower than
    `jacobian_perturbation_tol`, we increase it until its absolute value is higher than
    `jacobian_perturbation_tol`.
    (**Default**: 1e-7)
- `max_iterations::Int`: Maximum number of iterations allowed for the least-square fitting.
    (**Default**: 50)
- `mean_elements_epoch::Number`: Epoch of the fitted mean elements [Julian Day, UTC].
    (**Default**: `vjd[end]`)
- `template::Union{Nothing, sink, NamedTuple}`: Source of the metadata of the output. If it
    is an object of type `sink`, its metadata is copied. If it is a `NamedTuple`, its
    entries are passed as keywords to the constructor of `sink` on top of the default
    metadata, e.g. `(; name = "AMAZONIA 1", satellite_number = 47699)` for a `TLE` or
    `(; object_name = "AMAZONIA 1", object_id = "2021-015A", norad_cat_id = 47699)` for
    an `OrbitMeanElementsMessage`.
    (**Default**: `nothing`)
- `verbose::Bool`: If `true`, the algorithm prints debugging information to `stdout`.
    (**Default**: `true`)
- `weight_vector::AbstractVector`: Vector with the measurements weights for the least-square
    algorithm. We assemble the weight matrix `W` as a diagonal matrix with the elements in
    `weight_vector` at its diagonal.
    (**Default**: `@SVector(ones(Bool, 6))`)

```@repl sgp4
vr_teme = [
    [-6792.402703741442, 2192.6458461287293, 0.18851758695295118] .* 1000,
    [-6357.88873265975, 2391.9476768911686, 2181.838771262736] .* 1000
];

vv_teme = [
    [0.3445760107690598, 1.0395135806993514, 7.393686131436984] .* 1000,
    [2.5285015912807003, 0.27812476784300005, 7.030323100703928] .* 1000
];

vjd = [
    2.46002818657856e6,
    2.460028190050782e6
];

omm, P = Propagators.fit_mean_elements(Val(:SGP4), vjd, vr_teme, vv_teme; estimate_bstar = false)

omm

tle, P = Propagators.fit_mean_elements(TLE, Val(:SGP4), vjd, vr_teme, vv_teme; estimate_bstar = false)

tle
```

## References

- **[1]** **Hoots, F. R., Roehrich, R. L** (1980). *Models for Propagation of NORAD Elements
  Set*. **Spacetrack Report No. 3**.
- **[2]** **Vallado, D. A., Crawford, P., Hujsak, R., Kelso, T. S** (2006). *Revisiting
  Spacetrack Report #3: Rev1*. **AIAA**.
- **[3]** SGP4 Source code of [STRF](https://github.com/cbassa/strf), which the C code was
  converted by Paul. S. Crawford and Andrew R. Brooks.
- **[4]** **Vallado, D. A., Crawford, P** (2008). *SGP4 Orbit Determination*. **AIAA**.
