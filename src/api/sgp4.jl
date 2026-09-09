## Description #############################################################################
#
#  API implementation for SGP4 orbit propagator.
#
############################################################################################

# Implement the `Propagators` API for the SGP4 orbit propagator.
Propagators.epoch(orbp::OrbitPropagatorSgp4)           = orbp.sgp4d.epoch
Propagators.is_initialized(orbp::OrbitPropagatorSgp4)  = _is_initialized(orbp.sgp4d)
Propagators.last_instant(orbp::OrbitPropagatorSgp4)    = orbp.sgp4d.Δt * 60
Propagators.name(orbp::OrbitPropagatorSgp4)            = "SGP4 Orbit Propagator"
Propagators.propagator_data(orbp::OrbitPropagatorSgp4) = orbp.sgp4d

"""
    Propagators.mean_elements(
        orbp::OrbitPropagatorSgp4{Tepoch, T}
    ) where {Tepoch <: Number, T <: Number} -> KeplerianElements{MeanAnomaly, Tepoch, T}

Return the mean Keplerian elements [SI units] of the SGP4 orbit propagator `orbp` at the
last propagation instant. The initial TLE is rebuilt from the propagator and its epoch is
updated to that instant, which fails for epochs outside the years a TLE can represent.

# Extended help

## Throws

- `ArgumentError`: If the epoch year is outside the interval [1976, 2075].
"""
function Propagators.mean_elements(orbp::OrbitPropagatorSgp4)
    # We need to copy the propagator because updating the TLE epoch modifies it.
    sgp4d = copy(orbp.sgp4d)
    sgp4c = sgp4d.sgp4c

    # First, we need to create a TLE based on the initial parameters.
    dt  = julian2datetime(Propagators.epoch(orbp))
    dt₀ = DateTime(Year(dt))

    # A TLE stores the epoch year with two digits, which are interpreted as 1900 + y if
    # y > 75 and as 2000 + y otherwise. Hence, only the years between 1976 and 2075 can be
    # represented.
    dt_year = year(dt)

    if (dt_year < 1976) || (dt_year > 2075)
        throw(
            ArgumentError(
                "The epoch year $dt_year cannot be represented in a TLE, which only supports " *
                "the years between 1976 and 2075.",
            ),
        )
    end

    epoch_year = mod(dt_year, 100)
    epoch_day  = (dt - dt₀).value / 1000 / 86400 + 1

    tle = TLE(;
        epoch_year          = epoch_year,
        epoch_day           = epoch_day,
        bstar               = sgp4d.bstar,
        inclination         = sgp4d.i₀ |> rad2deg,
        raan                = sgp4d.Ω₀ |> rad2deg,
        eccentricity        = sgp4d.e₀,
        argument_of_perigee = sgp4d.ω₀ |> rad2deg,
        mean_anomaly        = sgp4d.M₀ |> rad2deg,
        mean_motion         = 720 * sgp4d.n₀ / π,
    )

    # Now, we update the TLE epoch to the current propagation instant, which also
    # initializes the copied propagator with the updated mean elements.
    new_epoch = Propagators.epoch(orbp) + Propagators.last_instant(orbp) / 86400
    update_sgp4_mean_elements_epoch!(sgp4d, tle, new_epoch; verbose = false)

    # Create and return the Keplerian elements storing the mean anomaly. The semi-major axis
    # is recovered from the mean motion using the SGP4 constants, converting it to meters.
    return KeplerianElements{MeanAnomaly}(
        new_epoch,
        (sgp4c.XKE / sgp4d.n₀)^(2 // 3) * (1000 * sgp4c.R0),
        sgp4d.e₀,
        sgp4d.i₀,
        sgp4d.Ω₀,
        sgp4d.ω₀,
        sgp4d.M₀,
    )
end

"""
    Propagators.fit_mean_elements(
        [sink::Type, ]::Val{:SGP4},
        vjd::AbstractVector{Tjd},
        vr_teme::AbstractVector{Tv},
        vv_teme::AbstractVector{Tv};
        kwargs...
    ) where {Tjd <: Number, Tv <: AbstractVector} -> sink, SMatrix{7, 7, T}, NamedTuple

Fit a set of SGP4 mean elements using the osculating elements represented by a set of
position vectors `vr_teme` [m] and a set of velocity vectors `vv_teme` [m / s] represented
in the True-Equator, Mean-Equinox reference frame (TEME) at instants in the array `vjd`
[Julian Day, UTC]. The mean elements are returned as an object of type `sink`, which can be
`TLE` or `OrbitMeanElementsMessage`. If it is omitted, an `OrbitMeanElementsMessage` is
returned. The fitting can fail if the least-square iterations diverge.

This algorithm was based on **[1]**.

!!! note

    This algorithm version will allocate a new SGP4 propagator with the default constants
    `SGP4C_WGS84`. If another set of constants are required, use the function
    [`Propagators.fit_mean_elements!`](@ref) instead.

# Keywords

- `atol::Number`: Tolerance for the residual absolute value. If the residual is lower than
    `atol` at any iteration, the computation loop stops.
    (**Default**: 2e-4)
- `rtol::Number`: Tolerance for the relative difference between the residuals. If the
    relative difference between the residuals in two consecutive iterations is lower than
    `rtol`, the computation loop stops.
    (**Default**: 2e-4)
- `estimate_bstar::Bool`: If `true`, the algorithm will try to estimate the B* parameter.
    Otherwise, it will be set to 0 or to the value in the initial guess (see section
    **Initial Guess**).
    (**Default**: `true`)
- `include_covariance::Bool`: If `true`, the covariance of the mean position and velocity
    obtained by the least-square algorithm is stored in the covariance matrix section of
    the output, represented in the TEME reference frame. It is only used when `sink` is
    `OrbitMeanElementsMessage`.
    (**Default**: `true`)
- `initial_guess::Union{Nothing, AbstractVector, TLE, OrbitMeanElementsMessage}`: Initial
    guess for the fitting process. If it is `nothing`, the algorithm will obtain an initial
    estimate from the osculating elements in `vr_teme` and `vv_teme`. For more information,
    see the section **Initial Guess**.
    (**Default**: `nothing`)
- `jacobian_method::AbstractJacobianMethod`: Method used to compute the Jacobian matrix. It
    can be `FiniteDiffJacobian()` for finite differences or `ForwardDiffJacobian()` for
    **ForwardDiff.jl** automatic differentiation.
    (**Default**: `FiniteDiffJacobian()`)
- `jacobian_perturbation::Number`: Initial state perturbation to compute the
    finite-difference when calculating the Jacobian matrix. Only used with
    `FiniteDiffJacobian()`.
    (**Default**: 1e-3)
- `jacobian_perturbation_tol::Number`: Tolerance to accept the perturbation when calculating
    the Jacobian matrix. If the computed perturbation is lower than
    `jacobian_perturbation_tol`, we increase it until its absolute value is higher than
    `jacobian_perturbation_tol`. Only used with `FiniteDiffJacobian()`.
    (**Default**: 1e-7)
- `max_iterations::Int`: Maximum number of iterations allowed for the least-square fitting.
    (**Default**: 50)
- `mean_elements_epoch::Union{Number, DateTime}`: Epoch of the fitted mean elements,
    represented by a Julian Day [UTC] or a `DateTime` [UTC].
    (**Default**: `vjd[end]`)
- `template::Union{Nothing, sink, NamedTuple}`: Source of the metadata of the output. If it
    is an object of type `sink`, its metadata is copied, e.g. the satellite name and number
    of a `TLE` or the header, the metadata, and the TLE-related parameters of an
    `OrbitMeanElementsMessage`. If it is a `NamedTuple`, its entries are passed as keywords
    to the constructor of `sink` on top of the default metadata, e.g.
    `(; name = "AMAZONIA 1", satellite_number = 47699)` for a `TLE` or
    `(; object_name = "AMAZONIA 1", object_id = "2021-015A", norad_cat_id = 47699)` for
    an `OrbitMeanElementsMessage`. In both cases, only the mean elements, the epoch, the
    drag term, and the covariance matrix are set by the fitted values, and a `NamedTuple`
    must not contain them. If it is `nothing`, the metadata is filled with default values.
    (**Default**: `nothing`)
- `verbose::Bool`: If `true`, the algorithm prints debugging information to `stdout`.
    (**Default**: `true`)
- `weight_vector::AbstractVector`: Vector with the measurements weights for the least-square
    algorithm. We assemble the weight matrix `W` as a diagonal matrix with the elements in
    `weight_vector` at its diagonal.
    (**Default**: `@SVector(ones(Bool, 6))`)

# Returns

- `sink`: The fitted mean elements.
- `SMatrix{7, 7, T}`: Final covariance matrix of the least-square algorithm, whose state is
    the mean position [km], the mean velocity [km / s], and the drag term B* [1 / er].
- `NamedTuple`: Statistics of the least-square algorithm with the following fields:
    - `converged::Bool`: `true` if the iterations stopped because the residue was lower
        than `atol` or its relative variation was lower than `rtol`, or `false` if they
        stopped by reaching `max_iterations`.
    - `iterations::Int`: Number of iterations performed.
    - `position_rmse::T`: RMSE of the position residue in the last iteration [m].
    - `velocity_rmse::T`: RMSE of the velocity residue in the last iteration [m / s].
    - `total_rmse::T`: Weighted RMSE of the residue in the last iteration, scaled to SI
        units.

    The statistics refer to the fitting of the mean elements. If their epoch is updated
    afterward to match `mean_elements_epoch`, the statistics of that update are not
    returned. The RMSE values are computed by **SatelliteToolboxSgp4.jl** in kilometers
    and converted here to SI units.

# Initial Guess

This algorithm uses a least-square algorithm to fit a set of mean elements based on a set
of osculating state vectors. Since the system is chaotic, a good initial guess is paramount
for algorithm convergence. We can provide an initial guess using the keyword
`initial_guess`.

If `initial_guess` is a `TLE` or an `OrbitMeanElementsMessage`, we update its epoch to the
desired one in `mean_elements_epoch`. Afterward, we use these mean elements as the initial
guess.

If `initial_guess` is an `AbstractVector`, we use this vector as the initial mean state
vector for the algorithm. It must contain 7 elements as follows:

    ┌                                    ┐
    │ IDs 1 to 3: Mean position [km]     │
    │ IDs 4 to 6: Mean velocity [km / s] │
    │ ID  7:      Bstar         [1 / er] │
    └                                    ┘

If `initial_guess` is `nothing`, the algorithm takes the closest osculating state vector to
the `mean_elements_epoch` and uses it as the initial mean state vector. In this case, the
epoch is set to the same epoch of the osculating data in `vjd`. When the fitted mean
elements are obtained, the algorithm updates their epoch to `mean_elements_epoch`.

!!! note

    If `initial_guess` is not `nothing`, the B* initial estimate is obtained from the mean
    elements or the state vector. Hence, if `estimate_bstar` is `false`, it will be kept
    constant with this initial value.

# References

- **[1]** Vallado, D. A., Crawford, P (2008). SGP4 Orbit Determination. American Institute
    of Aeronautics and Astronautics.

# Extended help

## Throws

- `ArgumentError`: If `vjd`, `vr_teme`, and `vv_teme` do not have the same length, if
    `weight_vector` does not have six elements, if `max_iterations` is lower than 1, or if
    a `NamedTuple` template contains a field set by the fit.
- `Sgp4FitDivergenceError`: If the least-square iterations diverge.
"""
function Propagators.fit_mean_elements(
    ::Val{:SGP4},
    vjd::AbstractVector{Tjd},
    vr_teme::AbstractVector{Tv},
    vv_teme::AbstractVector{Tv};
    kwargs...,
) where {Tjd <: Number, Tv <: AbstractVector}
    return Propagators.fit_mean_elements(
        OrbitMeanElementsMessage, Val(:SGP4), vjd, vr_teme, vv_teme; kwargs...
    )
end

function Propagators.fit_mean_elements(
    ::Type{S},
    ::Val{:SGP4},
    vjd::AbstractVector{Tjd},
    vr_teme::AbstractVector{Tv},
    vv_teme::AbstractVector{Tv};
    kwargs...,
) where {S <: Union{TLE, OrbitMeanElementsMessage}, Tjd <: Number, Tv <: AbstractVector}
    me, P, stats = fit_sgp4_mean_elements(
        S, vjd, vr_teme ./ 1000, vv_teme ./ 1000; kwargs...
    )
    return me, P, _stats_to_si(stats)
end

"""
    Propagators.fit_mean_elements!(
        orbp::OrbitPropagatorSgp4{Tepoch, T},
        vjd::AbstractVector{Tjd},
        vr_teme::AbstractVector{Tv},
        vv_teme::AbstractVector{Tv}[, sink::Type];
        kwargs...
    ) where {
        Tepoch <: Number,
        T <: Number,
        Tjd <: Number,
        Tv <: AbstractVector
    } -> sink, SMatrix{7, 7, T}, NamedTuple

Fit a set of SGP4 mean elements for the orbit propagator `orbp` using the osculating
elements represented by a set of position vectors `vr_teme` [m] and a set of velocity
vectors `vv_teme` [m / s] represented in the True-Equator, Mean-Equinox reference frame
(TEME) at instants in the array `vjd` [Julian Day, UTC]. The mean elements are returned as
an object of type `sink`, which can be `TLE` or `OrbitMeanElementsMessage`. If it is
omitted, an `OrbitMeanElementsMessage` is returned. The fitting can fail if the
least-square iterations diverge.

This algorithm was based on **[1]**.

!!! note

    The SGP4 orbit propagator `orbp` will be initialized with the mean elements returned by
    the function.

# Keywords

The keywords are the same as in [`Propagators.fit_mean_elements`](@ref) for the SGP4 orbit
propagator.

# Returns

- `sink`: The fitted mean elements.
- `SMatrix{7, 7, T}`: Final covariance matrix of the least-square algorithm, whose state is
    the mean position [km], the mean velocity [km / s], and the drag term B* [1 / er].
- `NamedTuple`: Statistics of the least-square algorithm with the following fields:
    - `converged::Bool`: `true` if the iterations stopped because the residue was lower
        than `atol` or its relative variation was lower than `rtol`, or `false` if they
        stopped by reaching `max_iterations`.
    - `iterations::Int`: Number of iterations performed.
    - `position_rmse::T`: RMSE of the position residue in the last iteration [m].
    - `velocity_rmse::T`: RMSE of the velocity residue in the last iteration [m / s].
    - `total_rmse::T`: Weighted RMSE of the residue in the last iteration, scaled to SI
        units.

    The statistics refer to the fitting of the mean elements. If their epoch is updated
    afterward to match `mean_elements_epoch`, the statistics of that update are not
    returned. The RMSE values are computed by **SatelliteToolboxSgp4.jl** in kilometers
    and converted here to SI units.

# References

- **[1]** Vallado, D. A., Crawford, P (2008). SGP4 Orbit Determination. American Institute
    of Aeronautics and Astronautics.

# Extended help

## Throws

- `ArgumentError`: If `vjd`, `vr_teme`, and `vv_teme` do not have the same length, if
    `weight_vector` does not have six elements, if `max_iterations` is lower than 1, or if
    a `NamedTuple` template contains a field set by the fit.
- `Sgp4FitDivergenceError`: If the least-square iterations diverge.
"""
function Propagators.fit_mean_elements!(
    orbp::OrbitPropagatorSgp4,
    vjd::AbstractVector{Tjd},
    vr_teme::AbstractVector{Tv},
    vv_teme::AbstractVector{Tv};
    kwargs...,
) where {Tjd <: Number, Tv <: AbstractVector}
    return Propagators.fit_mean_elements!(
        orbp, vjd, vr_teme, vv_teme, OrbitMeanElementsMessage; kwargs...
    )
end

function Propagators.fit_mean_elements!(
    orbp::OrbitPropagatorSgp4,
    vjd::AbstractVector{Tjd},
    vr_teme::AbstractVector{Tv},
    vv_teme::AbstractVector{Tv},
    ::Type{S};
    kwargs...,
) where {S <: Union{TLE, OrbitMeanElementsMessage}, Tjd <: Number, Tv <: AbstractVector}
    me, P, stats = fit_sgp4_mean_elements!(
        orbp.sgp4d, S, vjd, vr_teme ./ 1000, vv_teme ./ 1000; kwargs...
    )
    return me, P, _stats_to_si(stats)
end

"""
    Propagators.init(
        Val(:SGP4),
        epoch::Number,
        n₀::Number,
        e₀::Number,
        i₀::Number,
        Ω₀::Number,
        ω₀::Number,
        M₀::Number,
        bstar::Number;
        kwargs...
    ) -> OrbitPropagatorSgp4

    Propagators.init(Val(:SGP4), tle::TLE; kwargs...) -> OrbitPropagatorSgp4

    Propagators.init(
        Val(:SGP4),
        omm::OrbitMeanElementsMessage;
        kwargs...
    ) -> OrbitPropagatorSgp4

Create and initialize the SGP4 orbit propagator structure using the initial orbit specified
by the arguments. The initialization from an Orbit Mean-Elements Message (OMM) can fail if
the message does not describe an SGP4 orbit.

!!! note

    The type used in the propagation will be the same as used to define the constants in the
    structure `sgp4c`.

# Arguments

- `epoch::Number`: Epoch of the orbital elements [Julian Day].
- `n₀::Number`: SGP type "mean" mean motion at epoch [rad/s].
- `e₀::Number`: "Mean" eccentricity at epoch.
- `i₀::Number`: "Mean" inclination at epoch [rad].
- `Ω₀::Number`: "Mean" longitude of the ascending node at epoch [rad].
- `ω₀::Number`: "Mean" argument of perigee at epoch [rad].
- `M₀::Number`: "Mean" mean anomaly at epoch [rad].
- `bstar::Number`: Drag parameter (B*).
- `tle::TLE`: Two-line elements used for the initialization.
- `omm::OrbitMeanElementsMessage`: Orbit Mean-Elements Message used for the initialization.
    Its mean element theory must be `"SGP4"`, and it must provide the mean motion or the
    semi-major axis together with the gravitational coefficient. The drag term is set to 0
    if it is absent. See `sgp4_init` in **SatelliteToolboxSgp4.jl** for the details.

# Keywords

- `sgp4c::Sgp4Constants`: SGP4 orbit propagator constants (see `Sgp4Constants`). The
    constants in another number type can be obtained with the converting constructor, e.g.
    `Sgp4Constants{Float32}(SGP4C_WGS84)`.
    (**Default**: `SGP4C_WGS84`)

# Extended help

## Throws

- `ArgumentError`: If the mean element theory of `omm` is not `"SGP4"`, or if `omm`
    provides neither the mean motion nor the semi-major axis together with the
    gravitational coefficient.
"""
function Propagators.init(::Val{:SGP4}, tle::TLE; sgp4c::Sgp4Constants = SGP4C_WGS84)
    sgp4d = sgp4_init(tle; sgp4c = sgp4c)
    return OrbitPropagatorSgp4(sgp4d)
end

function Propagators.init(
    ::Val{:SGP4}, omm::OrbitMeanElementsMessage; sgp4c::Sgp4Constants = SGP4C_WGS84
)
    sgp4d = sgp4_init(omm; sgp4c = sgp4c)
    return OrbitPropagatorSgp4(sgp4d)
end

function Propagators.init(
    ::Val{:SGP4},
    epoch::Number,
    n₀::Number,
    e₀::Number,
    i₀::Number,
    Ω₀::Number,
    ω₀::Number,
    M₀::Number,
    bstar::Number;
    sgp4c::Sgp4Constants = SGP4C_WGS84,
)
    sgp4d = sgp4_init(epoch, 60n₀, e₀, i₀, Ω₀, ω₀, M₀, bstar; sgp4c = sgp4c)
    return OrbitPropagatorSgp4(sgp4d)
end

"""
    Propagators.init!(
        orbp::OrbitPropagatorSgp4,
        epoch::Number,
        n₀::Number,
        e₀::Number,
        i₀::Number,
        Ω₀::Number,
        ω₀::Number,
        M₀::Number,
        bstar::Number;
        kwargs...
    ) -> Nothing

    Propagators.init!(
        orbp::OrbitPropagatorSgp4,
        tle::TLE;
        kwargs...
    ) -> Nothing

    Propagators.init!(
        orbp::OrbitPropagatorSgp4,
        omm::OrbitMeanElementsMessage;
        kwargs...
    ) -> Nothing

Initialize the SGP4 orbit propagator structure `orbp` using the initial orbit specified
by the arguments. The initialization from an Orbit Mean-Elements Message (OMM) can fail if
the message does not describe an SGP4 orbit.

!!! warning

    The propagation constants `sgp4c::Sgp4Constants` in `orbp.sgp4d` will not be changed.
    Hence, they must be initialized, e.g. by creating the structure with
    `Sgp4Propagator(SGP4C_WGS84)`.

# Arguments

- `epoch::Number`: Epoch of the orbital elements [Julian Day].
- `n₀::Number`: SGP type "mean" mean motion at epoch [rad/s].
- `e₀::Number`: "Mean" eccentricity at epoch.
- `i₀::Number`: "Mean" inclination at epoch [rad].
- `Ω₀::Number`: "Mean" longitude of the ascending node at epoch [rad].
- `ω₀::Number`: "Mean" argument of perigee at epoch [rad].
- `M₀::Number`: "Mean" mean anomaly at epoch [rad].
- `bstar::Number`: Drag parameter (B*).
- `tle::TLE`: Two-line elements used for the initialization.
- `omm::OrbitMeanElementsMessage`: Orbit Mean-Elements Message used for the initialization.
    Its mean element theory must be `"SGP4"`, and it must provide the mean motion or the
    semi-major axis together with the gravitational coefficient. The drag term is set to 0
    if it is absent. See `sgp4_init!` in **SatelliteToolboxSgp4.jl** for the details.

# Extended help

## Throws

- `ArgumentError`: If the mean element theory of `omm` is not `"SGP4"`, or if `omm`
    provides neither the mean motion nor the semi-major axis together with the
    gravitational coefficient.
"""
function Propagators.init!(orbp::OrbitPropagatorSgp4, tle::TLE)
    sgp4_init!(orbp.sgp4d, tle)
    return nothing
end

function Propagators.init!(orbp::OrbitPropagatorSgp4, omm::OrbitMeanElementsMessage)
    sgp4_init!(orbp.sgp4d, omm)
    return nothing
end

function Propagators.init!(
    orbp::OrbitPropagatorSgp4,
    epoch::Number,
    n₀::Number,
    e₀::Number,
    i₀::Number,
    Ω₀::Number,
    ω₀::Number,
    M₀::Number,
    bstar::Number,
)
    sgp4_init!(orbp.sgp4d, epoch, 60n₀, e₀, i₀, Ω₀, ω₀, M₀, bstar)
    return nothing
end

"""
    Propagators.propagate!(
        orbp::OrbitPropagatorSgp4{Tepoch, T},
        Δt::Number
    ) where {Tepoch <: Number, T <: Number} -> SVector{3, T}, SVector{3, T}

Propagate the orbit of the SGP4 orbit propagator `orbp` to `Δt` [s] after the epoch of the
TLE, updating the internal state of `orbp`. The SGP4 kernel works in minutes and kilometers,
so the instant and the output are converted to SI units here.

# Returns

- `SVector{3, T}`: Position vector [m] represented in the TEME frame at propagation instant.
- `SVector{3, T}`: Velocity vector [m / s] represented in the TEME frame at propagation
    instant.
"""
function Propagators.propagate!(orbp::OrbitPropagatorSgp4, Δt::Number)
    r_teme, v_teme = sgp4!(orbp.sgp4d, Δt / 60)
    return 1000r_teme, 1000v_teme
end

############################################################################################
#                                        Julia API                                         #
############################################################################################

function Base.copy(
    orbp::OrbitPropagatorSgp4{Tepoch, T}
) where {Tepoch <: Number, T <: Number}
    return OrbitPropagatorSgp4{Tepoch, T}(copy(orbp.sgp4d))
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

"""
    _stats_to_si(stats::NamedTuple) -> NamedTuple

Convert the statistics `stats` returned by the SGP4 fitting functions of
**SatelliteToolboxSgp4.jl**, whose RMSE values are in kilometers and kilometers per second,
to SI units.
"""
function _stats_to_si(stats::NamedTuple)
    return (;
        stats.converged,
        stats.iterations,
        position_rmse = 1000 * stats.position_rmse,
        velocity_rmse = 1000 * stats.velocity_rmse,
        total_rmse    = 1000 * stats.total_rmse,
    )
end
