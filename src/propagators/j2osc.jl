## Description #############################################################################
#
# J2 osculating orbit propagator algorithm.
#
# This algorithm propagates the orbit considering the secular and short-period
# perturbations introduced by the J2 gravitational term. The algorithm is based on Kwok
# version as indicated in [1, p. 708-710].
#
## References ##############################################################################
#
# [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. Microcosm
#     Press, Hawthorn, CA, USA.
#
############################################################################################

export j2osc_init, j2osc_init!, j2osc, j2osc!
export fit_j2osc_mean_elements, fit_j2osc_mean_elements!
export update_j2osc_mean_elements_epoch, update_j2osc_mean_elements_epoch!

############################################################################################
#                                        Julia API                                         #
############################################################################################

# Define `copy` for the propagator structure.
let
    fields = fieldnames(J2OsculatingPropagator)
    expressions = vcat(
        :(new_j2oscd.j2d = copy(j2oscd.j2d)),
        [:(new_j2oscd.$f = j2oscd.$f) for f in fields if f != :j2d],
    )

    @eval begin
        function Base.copy(
            j2oscd::J2OsculatingPropagator{Tepoch, T}
        ) where {Tepoch <: Number, T <: Number}
            new_j2oscd = J2OsculatingPropagator{Tepoch, T}()
            $(expressions...)
            return new_j2oscd
        end
    end
end

############################################################################################
#                                        Functions                                         #
############################################################################################

"""
    j2osc_init(orb₀::KeplerianElements; kwargs...) -> J2OsculatingPropagator

Create and initialize the J2 osculating orbit propagator structure using the mean Keplerian
elements `orb₀` [SI units].

!!! note

    The type used in the propagation will be the same as used to define the constants in the
    structure `j2c`.

# Keywords

- `j2c::J2PropagatorConstants`: J2 orbit propagator constants (see
    [`J2PropagatorConstants`](@ref)).
    (**Default** = `j2c_egm2008`)
"""
function j2osc_init(
    orb₀::KeplerianElements{Tepoch, Tkepler}; j2c::J2PropagatorConstants{T} = j2c_egm2008
) where {Tepoch <: Number, Tkepler <: AbstractFloat, T <: Number}
    # Allocate the J2 propagator structure that will propagate the mean elements.
    j2d = J2Propagator{Tepoch, T}()

    # Assign the constants, which are used in the initialization.
    j2d.j2c = j2c

    # Allocate the J2 osculating propagator structure.
    j2oscd = J2OsculatingPropagator{Tepoch, T}()
    j2oscd.j2d = j2d

    # Initialize the propagator and return.
    j2osc_init!(j2oscd, orb₀)

    return j2oscd
end

function j2osc_init(
    orb₀::KeplerianElements{Tepoch, Tkepler}; j2c::J2PropagatorConstants{Tj2c} = j2c_egm2008
) where {Tepoch <: Number, Tkepler <: Number, Tj2c <: Number}
    T = promote_type(Tj2c, Tkepler)

    # Allocate the J2 propagator structure that will propagate the mean elements.
    j2d = J2Propagator{Tepoch, T}()

    # Assign the constants, which are used in the initialization.
    j2d.j2c = J2PropagatorConstants{T}(j2c.R0, j2c.μm, j2c.J2)

    # Allocate the J2 osculating propagator structure.
    j2oscd = J2OsculatingPropagator{Tepoch, T}()
    j2oscd.j2d = j2d

    # Initialize the propagator and return.
    j2osc_init!(j2oscd, orb₀)

    return j2oscd
end

"""
    j2osc_init!(j2oscd::J2OsculatingPropagator, orb₀::KeplerianElements) -> Nothing

Initialize the J2 osculating orbit propagator structure `j2oscd` using the mean Keplerian
elements `orb₀` [SI units].

!!! warning

    The propagation constants `j2c::J2PropagatorConstants` in `j2oscd.j2d` will not be
    changed. Hence, they must be initialized.
"""
function j2osc_init!(j2oscd::J2OsculatingPropagator, orb₀::KeplerianElements)
    _j2osc_init!(j2oscd, orb₀)

    # Call the propagation one time to update the osculating elements.
    j2osc!(j2oscd, 0)

    return nothing
end

# Initialize the propagator without computing the osculating elements at the initial instant.
# The callers that propagate the orbit right afterwards use this function, since the
# osculating elements would be overwritten anyway.
function _j2osc_init!(j2oscd::J2OsculatingPropagator, orb₀::KeplerianElements)
    # Initialize the J2 propagator that will propagate the mean elements.
    j2_init!(j2oscd.j2d, orb₀)

    return nothing
end

"""
    j2osc(Δt::Number, orb₀::KeplerianElements; kwargs...) -> SVector{3, T}, SVector{3, T}, J2OsculatingPropagator

Initialize the J2 osculating propagator structure using the input elements `orb₀` [SI units]
and propagate the orbit until the time Δt [s].

!!! note

    The type used in the propagation will be the same as used to define the constants in the
    structure `j2c`.

# Keywords

- `j2c::J2PropagatorConstants{T}`: J2 orbit propagator constants (see
    [`J2PropagatorConstants`](@ref)).
    (**Default** = `j2c_egm2008`)

# Returns

- `SVector{3, T}`: Position vector [m] represented in the inertial frame at propagation
    instant.
- `SVector{3, T}`: Velocity vector [m / s] represented in the inertial frame at propagation
    instant.
- [`J2OsculatingPropagator`](@ref): Structure with the initialized propagator.

# Remarks

The inertial frame in which the output is represented depends on which frame was used to
generate the orbit parameters. Notice that the perturbation theory requires an inertial
frame with true equator.
"""
function j2osc(
    Δt::Number, orb₀::KeplerianElements; j2c::J2PropagatorConstants = j2c_egm2008
)
    j2oscd = j2osc_init(orb₀; j2c = j2c)
    r_i, v_i = j2osc!(j2oscd, Δt)
    return r_i, v_i, j2oscd
end

"""
    j2osc!(j2oscd::J2OsculatingPropagator{Tepoch, T}, t::Number) where {Tepoch, T} -> SVector{3, T}, SVector{3, T}

Propagate the orbit defined in `j2oscd` (see [`J2OsculatingPropagator`](@ref)) to `t` [s]
after the epoch of the input mean elements in `j2d`.

!!! note

    The internal values in `j2oscd` will be modified.

# Returns

- `SVector{3, T}`: Position vector [m] represented in the inertial frame at propagation
    instant.
- `SVector{3, T}`: Velocity vector [m / s] represented in the inertial frame at propagation
    instant.

# Remarks

The inertial frame in which the output is represented depends on which frame was used to
generate the orbit parameters. Notice that the perturbation theory requires an inertial
frame with true equator.
"""
function j2osc!(
    j2oscd::J2OsculatingPropagator{Tepoch, T}, t::Number
) where {Tepoch <: Number, T <: Number}
    # First, we need to propagate the mean elements since they are necessary to compute the
    # short-periodic perturbations.
    j2d = j2oscd.j2d
    mean_orbk = _j2_mean_elements!(j2d, t)

    # Unpack the propagator constants.
    j2c = j2d.j2c
    R₀  = j2c.R0
    μm  = j2c.μm
    J₂  = j2c.J2

    orbk = _osculating_elements(
        mean_orbk, j2d.M_k, j2d.orb₀.t + Tepoch(t) / 86400, R₀, μm, J₂
    )

    # Compute the position and velocity considering the osculating elements.
    r_i_k, v_i_k = kepler_to_rv(orbk)

    # Update the J2 orbit propagator structure.
    j2oscd.Δt   = T(t)
    j2oscd.orbk = orbk

    return r_i_k, v_i_k
end

"""
    fit_j2osc_mean_elements(vjd::AbstractVector{Tjd}, vr_i::AbstractVector{Tv}, vv_i::AbstractVector{Tv}; kwargs...) where {Tjd<:Number, Tv<:AbstractVector} -> KeplerianElements{Float64, Float64}, SMatrix{6, 6, Float64}

Fit a set of mean Keplerian elements for the J2 osculating orbit propagator using the
osculating elements represented by a set of position vectors `vr_i` [m] and a set of
velocity vectors `vv_i` [m / s] represented in an inertial reference frame at instants in
the array `vjd` [Julian Day].

!!! note

    This algorithm version will allocate a new J2 osculating propagator with the default
    constants `j2c_egm2008`. If another set of constants are required, use the function
    [`fit_j2osc_mean_elements!`](@ref) instead.

# Keywords

- `atol::Number`: Tolerance for the residual absolute value. If the residual is lower than
    `atol` at any iteration, the computation loop stops.
    (**Default** = 2e-4)
- `rtol::Number`: Tolerance for the relative difference between the residuals. If the
    relative difference between the residuals in two consecutive iterations is lower than
    `rtol`, the computation loop stops.
    (**Default** = 2e-4)
- `initial_guess::Union{Nothing, KeplerianElements}`: Initial guess for the mean elements
    fitting process. If it is `nothing`, the algorithm will obtain an initial estimate from
    the osculating elements in `vr_i` and `vv_i`.
    (**Default** = nothing)
- `jacobian_method::Union{FiniteDiffJacobian, ForwardDiffJacobian}`: Method used to compute
    the Jacobian matrix. Use `FiniteDiffJacobian()` for finite differences or
    `ForwardDiffJacobian()` for `ForwardDiff.jl` automatic differentiation.
    (**Default** = `FiniteDiffJacobian()`)
- `jacobian_perturbation::Number`: Initial state perturbation to compute the
    finite-difference when calculating the Jacobian matrix. Only used with `FiniteDiffJacobian()`.
    (**Default** = 1e-3)
- `jacobian_perturbation_tol::Number`: Tolerance to accept the perturbation when calculating
    the Jacobian matrix. If the computed perturbation is lower than
    `jacobian_perturbation_tol`, we increase it until its absolute value is higher than
    `jacobian_perturbation_tol`. Only used with `FiniteDiffJacobian()`.
    (**Default** = 1e-7)
- `max_iterations::Int`: Maximum number of iterations allowed for the least-square fitting.
    (**Default** = 50)
- `mean_elements_epoch::Number`: Epoch for the fitted mean elements.
    (**Default** = vjd[end])
- `verbose::Bool`: If `true`, the algorithm prints debugging information to `stdout`.
    (**Default** = true)
- `weight_vector::AbstractVector`: Vector with the measurements weights for the least-square
    algorithm. We assemble the weight matrix `W` as a diagonal matrix with the elements in
    `weight_vector` at its diagonal.
    (**Default** = `@SVector(ones(Bool, 6))`)

# Returns

- `KeplerianElements{Float64, Float64}`: Fitted Keplerian elements.
- `SMatrix{6, 6, Float64}`: Final covariance matrix of the least-square algorithm.

# Examples

```julia-repl
julia> vr_i = [
           [-6792.402703741442, 2192.6458461287293, 0.18851758695295118] .* 1000,
           [-6357.88873265975, 2391.9476768911686, 2181.838771262736] .* 1000
       ];

julia> vv_i = [
           [0.3445760107690598, 1.0395135806993514, 7.393686131436984] .* 1000,
           [2.5285015912807003, 0.27812476784300005, 7.030323100703928] .* 1000
       ];

julia> vjd = [
           2.46002818657856e6,
           2.460028190050782e6
       ];

julia> orb, P = fit_j2osc_mean_elements(vjd, vr_i, vv_i)
ACTION:   Fitting the mean elements for the J2 osculating propagator.
           Iteration        Position RMSE        Velocity RMSE           Total RMSE       RMSE Variation
                                     [km]             [km / s]                  [ ]
PROGRESS:          4          1.69161e-05           0.00260193              2.60198         -2.06772e-09 %

(KeplerianElements{Float64, Float64}: Epoch = 2.46003e6 (2023-03-24T16:33:40.388), [0.9999427882949705 -0.0049011518111948425 … -0.00012021282848110985 -0.00011550570919063845; -0.00490115181520327 1.0007325785722665 … 0.0032661568999667063 3.8567115793316056e-5; … ; -0.00012021282849269182 0.003266156899966125 … 2.204826778451109e-5 5.382733068487529e-8; -0.00011550570919086411 3.8567115786674836e-5 … 5.382733066340212e-8 2.1657420198753813e-5])

julia> orb
KeplerianElements{Float64, Float64}:
           Epoch :    2.46003e6 (2023-03-24T16:33:40.388)
 Semi-major axis : 7135.8        km
    Eccentricity :    0.00135383
     Inclination :   98.4304     °
            RAAN :  162.113      °
 Arg. of Perigee :   64.9256     °
    True Anomaly :  313.085      °
```
"""
function fit_j2osc_mean_elements(
    vjd::AbstractVector{Tjd}, vr_i::AbstractVector{Tv}, vv_i::AbstractVector{Tv}; kwargs...
) where {Tjd <: Number, Tv <: AbstractVector}
    # Allocate the J2 propagator structure that will propagate the mean elements.
    j2d = J2Propagator{Float64, Float64}()

    # Assign the constants, which are used in the initialization.
    j2d.j2c = j2c_egm2008

    # Allocate the J2 osculating propagator structure.
    j2oscd = J2OsculatingPropagator{Float64, Float64}()
    j2oscd.j2d = j2d

    return fit_j2osc_mean_elements!(j2oscd, vjd, vr_i, vv_i; kwargs...)
end

"""
    fit_j2osc_mean_elements!(j2oscd::J2OsculatingPropagator{Tepoch, T}, vjd::AbstractVector{Tjd}, vr_i::AbstractVector{Tv}, vv_i::AbstractVector{Tv}; kwargs...) where {T<:Number, Tepoch<:Number, Tjd<:Number, Tv<:AbstractVector} -> KeplerianElements{Tepoch, T}, SMatrix{6, 6, T}

Fit a set of mean Keplerian elements for the J2 osculating orbit propagator `j2oscd` using
the osculating elements represented by a set of position vectors `vr_i` [m] and a set of
velocity vectors `vv_i` [m / s] represented in an inertial reference frame at instants in
the array `vjd` [Julian Day].

!!! note

    The J2 osculating orbit propagator `j2oscd` will be initialized with the Keplerian
    elements returned by the function.

# Keywords

- `atol::Number`: Tolerance for the residual absolute value. If the residual is lower than
    `atol` at any iteration, the computation loop stops.
    (**Default** = 2e-4)
- `rtol::Number`: Tolerance for the relative difference between the residuals. If the
    relative difference between the residuals in two consecutive iterations is lower than
    `rtol`, the computation loop stops.
    (**Default** = 2e-4)
- `initial_guess::Union{Nothing, KeplerianElements}`: Initial guess for the mean elements
    fitting process. If it is `nothing`, the algorithm will obtain an initial estimate from
    the osculating elements in `vr_i` and `vv_i`.
    (**Default** = nothing)
- `jacobian_method::Union{FiniteDiffJacobian, ForwardDiffJacobian}`: Method used to compute
    the Jacobian matrix. Use `FiniteDiffJacobian()` for finite differences or
    `ForwardDiffJacobian()` for `ForwardDiff.jl` automatic differentiation.
    (**Default** = `FiniteDiffJacobian()`)
- `jacobian_perturbation::Number`: Initial state perturbation to compute the
    finite-difference when calculating the Jacobian matrix. Only used with `FiniteDiffJacobian()`.
    (**Default** = 1e-3)
- `jacobian_perturbation_tol::Number`: Tolerance to accept the perturbation when calculating
    the Jacobian matrix. If the computed perturbation is lower than
    `jacobian_perturbation_tol`, we increase it until its absolute value is higher than
    `jacobian_perturbation_tol`. Only used with `FiniteDiffJacobian()`.
    (**Default** = 1e-7)
- `max_iterations::Int`: Maximum number of iterations allowed for the least-square fitting.
    (**Default** = 50)
- `mean_elements_epoch::Number`: Epoch for the fitted mean elements.
    (**Default** = vjd[end])
- `verbose::Bool`: If `true`, the algorithm prints debugging information to `stdout`.
    (**Default** = true)
- `weight_vector::AbstractVector`: Vector with the measurements weights for the least-square
    algorithm. We assemble the weight matrix `W` as a diagonal matrix with the elements in
    `weight_vector` at its diagonal.
    (**Default** = `@SVector(ones(Bool, 6))`)

# Returns

- `KeplerianElements{Tepoch, T}`: Fitted Keplerian elements.
- `SMatrix{6, 6, T}`: Final covariance matrix of the least-square algorithm.

# Examples

```julia-repl
# Allocate a new J2 osculating orbit propagator using a dummy set of Keplerian elements.
julia> j2oscd = j2osc_init(KeplerianElements{Float64, Float64}(0, 7000e3, 0, 0, 0, 0, 0));

julia> vr_i = [
           [-6792.402703741442, 2192.6458461287293, 0.18851758695295118] .* 1000,
           [-6357.88873265975, 2391.9476768911686, 2181.838771262736] .* 1000
       ];

julia> vv_i = [
           [0.3445760107690598, 1.0395135806993514, 7.393686131436984] .* 1000,
           [2.5285015912807003, 0.27812476784300005, 7.030323100703928] .* 1000
       ];

julia> vjd = [
           2.46002818657856e6,
           2.460028190050782e6
       ];

julia> orb, P = fit_j2osc_mean_elements!(j2oscd, vjd, vr_i, vv_i)
ACTION:   Fitting the mean elements for the J2 osculating propagator.
           Iteration        Position RMSE        Velocity RMSE           Total RMSE       RMSE Variation
                                     [km]             [km / s]                  [ ]
PROGRESS:          4          1.69161e-05           0.00260193              2.60198         -2.06772e-09 %

(KeplerianElements{Float64, Float64}: Epoch = 2.46003e6 (2023-03-24T16:33:40.388), [0.9999427882949705 -0.0049011518111948425 … -0.00012021282848110985 -0.00011550570919063845; -0.00490115181520327 1.0007325785722665 … 0.0032661568999667063 3.8567115793316056e-5; … ; -0.00012021282849269182 0.003266156899966125 … 2.204826778451109e-5 5.382733068487529e-8; -0.00011550570919086411 3.8567115786674836e-5 … 5.382733066340212e-8 2.1657420198753813e-5])

julia> orb
KeplerianElements{Float64, Float64}:
           Epoch :    2.46003e6 (2023-03-24T16:33:40.388)
 Semi-major axis : 7135.8        km
    Eccentricity :    0.00135383
     Inclination :   98.4304     °
            RAAN :  162.113      °
 Arg. of Perigee :   64.9256     °
    True Anomaly :  313.085      °
```
"""
function fit_j2osc_mean_elements!(
    j2oscd::J2OsculatingPropagator{Tepoch, T},
    vjd::AbstractVector{Tjd},
    vr_i::AbstractVector{Tv},
    vv_i::AbstractVector{Tv};
    kwargs...,
) where {Tepoch <: Number, T <: Number, Tjd <: Number, Tv <: AbstractVector}
    return _fit_mean_elements!(j2oscd, vjd, vr_i, vv_i; kwargs...)
end

"""
    update_j2osc_mean_elements_epoch(orb::KeplerianElements, new_epoch::Union{Number, DateTime}) -> KeplerianElements

Update the epoch of the mean elements `orb` using a J2 osculating orbit propagator to
`new_epoch`, which can be represented by a Julian Day or a `DateTime`.

!!! note

    This algorithm version will allocate a new J2 osculating propagator with the default
    constants `j2c_egm2008`. If another set of constants are required, use the function
    [`update_j2osc_mean_elements_epoch!`](@ref) instead.

# Examples

```julia-repl
julia> orb = KeplerianElements(
           DateTime("2023-01-01") |> datetime2julian,
           7190.982e3,
           0.001111,
           98.405 |> deg2rad,
           90     |> deg2rad,
           200    |> deg2rad,
           45     |> deg2rad
       )
KeplerianElements{Float64, Float64}:
           Epoch :    2.45995e6 (2023-01-01T00:00:00)
 Semi-major axis : 7190.98     km
    Eccentricity :    0.001111
     Inclination :   98.405    °
            RAAN :   90.0      °
 Arg. of Perigee :  200.0      °
    True Anomaly :   45.0      °

julia> update_j2osc_mean_elements_epoch(orb, DateTime("2023-01-02"))
KeplerianElements{Float64, Float64}:
           Epoch :    2.45995e6 (2023-01-02T00:00:00)
 Semi-major axis : 7190.98     km
    Eccentricity :    0.001111
     Inclination :   98.405    °
            RAAN :   90.9565   °
 Arg. of Perigee :  197.078    °
    True Anomaly :  127.291    °
```
"""
function update_j2osc_mean_elements_epoch(
    orb::KeplerianElements{Tepoch, T}, new_epoch::Union{Number, DateTime}
) where {T <: Number, Tepoch <: Number}
    # Allocate the J2 propagator structure that will propagate the mean elements.
    j2d = J2Propagator{Tepoch, T}()

    # Assign the constants, which are used in the initialization.
    j2d.j2c = j2c_egm2008

    # Allocate the J2 osculating propagator structure.
    j2oscd = J2OsculatingPropagator{Tepoch, T}()
    j2oscd.j2d = j2d

    return update_j2osc_mean_elements_epoch!(j2oscd, orb, new_epoch)
end

"""
    update_j2osc_mean_elements_epoch!(j2oscd::J2OsculatingPropagator, orb::KeplerianElements, new_epoch::Union{Number, DateTime}) -> KeplerianElements

Update the epoch of the mean elements `orb` using the propagator `j2oscd` to `new_epoch`,
which can be represented by a Julian Day or a `DateTime`.

!!! note

    The J2 osculating orbit propagator `j2oscd` will be initialized with the Keplerian
    elements returned by the function.

# Examples

```julia-repl
julia> orb = KeplerianElements(
           DateTime("2023-01-01") |> datetime2julian,
           7190.982e3,
           0.001111,
           98.405 |> deg2rad,
           90     |> deg2rad,
           200    |> deg2rad,
           45     |> deg2rad
       )
KeplerianElements{Float64, Float64}:
           Epoch :    2.45995e6 (2023-01-01T00:00:00)
 Semi-major axis : 7190.98     km
    Eccentricity :    0.001111
     Inclination :   98.405    °
            RAAN :   90.0      °
 Arg. of Perigee :  200.0      °
    True Anomaly :   45.0      °

# Allocate a new J2 osculating orbit propagator using the created Keplerian elements. Notice
# that any set of Keplerian elements can be used here.
julia> j2oscd = j2osc_init(orb);

julia> update_j2osc_mean_elements_epoch!(j2oscd, orb, DateTime("2023-01-02"))
KeplerianElements{Float64, Float64}:
           Epoch :    2.45995e6 (2023-01-02T00:00:00)
 Semi-major axis : 7190.98     km
    Eccentricity :    0.001111
     Inclination :   98.405    °
            RAAN :   90.9565   °
 Arg. of Perigee :  197.078    °
    True Anomaly :  127.291    °
```
"""
function update_j2osc_mean_elements_epoch!(
    j2oscd::J2OsculatingPropagator, orb::KeplerianElements, new_epoch::DateTime
)
    dt = datetime2julian(new_epoch)
    return update_j2osc_mean_elements_epoch!(j2oscd, orb, dt)
end

function update_j2osc_mean_elements_epoch!(
    j2oscd::J2OsculatingPropagator, orb::KeplerianElements, new_epoch::Number
)
    # First, we need to initialize the J2 osculating propagator with the mean elements.
    _j2osc_init!(j2oscd, orb)

    # Now, we just need to propagate the orbit to the desired instant and obtain the mean
    # elements from the J2 propagator structure inside.
    Δt = (new_epoch - j2oscd.j2d.orb₀.t) * 86400
    j2osc!(j2oscd, Δt)
    orb = j2oscd.j2d.orbk

    # Finally, we initialize the propagator with the new set of mean elements.
    j2osc_init!(j2oscd, orb)

    return orb
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

# Create a propagator that the finite-difference Jacobian can use as scratch space, so it
# does not clobber the propagator kept initialized by the fitting loop.
function _create_j2osc_fd_propagator(
    j2oscd::J2OsculatingPropagator{Tepoch, T}
) where {Tepoch, T}
    fd = J2OsculatingPropagator{Tepoch, T}()
    fd.j2d = J2Propagator{Tepoch, T}()
    fd.j2d.j2c = j2oscd.j2d.j2c
    return fd
end

function _create_j2osc_ad_propagator(
    j2oscd::J2OsculatingPropagator{Tepoch, T}
) where {Tepoch, T}
    tag = ForwardDiff.Tag{Nothing, T}
    D   = ForwardDiff.Dual{tag, T, 6}
    j2c = j2oscd.j2d.j2c

    ad = J2OsculatingPropagator{Tepoch, D}()
    ad.j2d = J2Propagator{Tepoch, D}()
    ad.j2d.j2c = J2PropagatorConstants{D}(D(j2c.R0), D(j2c.μm), D(j2c.J2))
    return ad
end

function _j2osc_jacobian(
    ::FiniteDiffJacobian,
    j2oscd::J2OsculatingPropagator{Tepoch, T},
    Δt::Number,
    x₁::SVector{6, T},
    y₁::SVector{6, T};
    perturbation::Number = T(1e-3),
    perturbation_tol::Number = T(1e-7),
    pd_ad::Union{Nothing, J2OsculatingPropagator} = nothing,
    pd_fd::Union{Nothing, J2OsculatingPropagator} = nothing,
) where {T <: Number, Tepoch <: Number}
    # The perturbed propagations below overwrite the propagator. Use the scratch propagator
    # when the caller provides one, so it can keep `j2oscd` initialized across the
    # measurements instead of reinitializing it for each one.
    epoch = j2oscd.j2d.orb₀.t
    fd::J2OsculatingPropagator{Tepoch, T} = isnothing(pd_fd) ? j2oscd : pd_fd

    J = MMatrix{6, 6, T}(undef)
    x₂ = x₁

    # Convert the perturbation parameters to the element type once. Otherwise, `ϵ` would be
    # assigned values of two different types, and the Jacobian columns would be computed in
    # the promoted type instead of in `T`.
    ϵ₀ = T(perturbation)
    ϵ_tol = T(perturbation_tol)

    @inbounds @views for j in 1:6
        α = x₂[j]
        ϵ = α * ϵ₀

        for _ in 1:5
            abs(ϵ) > ϵ_tol && break
            ϵ *= T(1.4)
        end

        if abs(ϵ) < ϵ_tol
            ϵ = signbit(α) ? -ϵ_tol : ϵ_tol
        end

        α += ϵ
        x₂ = setindex(x₂, α, j)

        orb = rv_to_kepler(x₂[SOneTo(3)], x₂[StaticArrays.SUnitRange(4, 6)], epoch)
        _j2osc_init!(fd, orb)
        r_i, v_i = j2osc!(fd, Δt)
        y₂ = @SVector [r_i[1], r_i[2], r_i[3], v_i[1], v_i[2], v_i[3]]

        J[:, j] .= (y₂ .- y₁) ./ ϵ
        x₂ = setindex(x₂, x₁[j], j)
    end

    return SMatrix{6, 6, T}(J)
end

function _j2osc_jacobian(
    ::ForwardDiffJacobian,
    j2oscd::J2OsculatingPropagator{Tepoch, T},
    Δt::Number,
    x₁::SVector{6, T},
    y₁::SVector{6, T};
    perturbation::Number = T(1e-3),
    perturbation_tol::Number = T(1e-7),
    pd_ad::Union{Nothing, J2OsculatingPropagator} = nothing,
    pd_fd::Union{Nothing, J2OsculatingPropagator} = nothing,
) where {T <: Number, Tepoch <: Number}
    epoch = j2oscd.j2d.orb₀.t
    N     = 6
    tag   = ForwardDiff.Tag{Nothing, T}
    D     = ForwardDiff.Dual{tag, T, N}

    # Declaring the type of the local makes the propagator concrete regardless of what
    # the caller provides. Otherwise, the unparameterised type in the keyword leaves
    # this call and the ones below dynamically dispatched.
    ad::J2OsculatingPropagator{Tepoch, D} =
        isnothing(pd_ad) ? _create_j2osc_ad_propagator(j2oscd) : pd_ad

    seeds  = ntuple(i -> ForwardDiff.Partials(ntuple(j -> T(i == j), Val(N))), Val(N))
    x_dual = SVector{N, D}(ntuple(i -> D(x₁[i], seeds[i]), Val(N)))

    orb = rv_to_kepler(x_dual[SOneTo(3)], x_dual[StaticArrays.SUnitRange(4, 6)], epoch)
    _j2osc_init!(ad, orb)
    r, v   = j2osc!(ad, Δt)
    y_dual = vcat(r, v)

    return SMatrix{6, N, T}(
        ntuple(k -> ForwardDiff.partials(y_dual[mod1(k, 6)], cld(k, 6)), Val(6 * N))
    )
end
