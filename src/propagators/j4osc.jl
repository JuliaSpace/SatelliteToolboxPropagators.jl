## Description #############################################################################
#
#   J4 osculating orbit propagator algorithm.
#
#   This algorithm propagates the orbit considering the secular perturbations considering
#   the terms J2, J2², and J4 and the short-period perturbations considering the J2
#   gravitation term only. The algorithm is based on Kwok version as indicated in
#   [1, p. 708-710].
#
## References ##############################################################################
#
# [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. 4th ed.
#     Microcosm Press, Hawthorn, CA, USA.
#
############################################################################################

export j4osc_init, j4osc_init!, j4osc, j4osc!
export fit_j4osc_mean_elements, fit_j4osc_mean_elements!
export update_j4osc_mean_elements_epoch, update_j4osc_mean_elements_epoch!

############################################################################################
#                                        Julia API                                         #
############################################################################################

# Define `copy` for the propagator structure.
let
    fields = fieldnames(J4OsculatingPropagator)
    expressions = vcat(
        :(new_j4oscd.j4d = copy(j4oscd.j4d)),
        [
            :(new_j4oscd.$f = j4oscd.$f)
            for f in fields if f != :j4d
        ]
    )

    @eval begin
        function Base.copy(
            j4oscd::J4OsculatingPropagator{Tepoch, T}
        ) where {Tepoch<:Number, T<:Number}
            new_j4oscd = J4OsculatingPropagator{Tepoch, T}()
            $(expressions...)
            return new_j4oscd
        end
    end
end

############################################################################################
#                                        Functions                                         #
############################################################################################

"""
    j4osc_init(orb₀::KeplerianElements; kwargs...) -> J4OsculatingPropagator

Create and initialize the J4 osculating orbit propagator structure using the mean Keplerian
elements `orb₀`.

!!! note

    The type used in the propagation will be the same as used to define the constants in the
    structure `j4c`.

# Keywords

- `j4c::J4PropagatorConstants`: J4 orbit propagator constants (see
    [`J4PropagatorConstants`](@ref)). (**Default** = `j4c_egm2008`)
"""
function j4osc_init(
    orb₀::KeplerianElements{Tepoch, Tkepler};
    j4c::J4PropagatorConstants{T} = j4c_egm2008
) where {Tepoch<:Number, Tkepler<:Number, T<:Number}
    # Allocate the J4 propagator structure that will propagate the mean elements.
    j4d = J4Propagator{Tepoch, T}()

    # Assign the constants, which are used in the initialization.
    j4d.j4c = j4c

    # Allocate the J4 osculating propagator structure.
    j4oscd = J4OsculatingPropagator{Tepoch, T}()
    j4oscd.j4d = j4d

    # Initialize the propagator and return.
    j4osc_init!(j4oscd, orb₀)

    return j4oscd
end

"""
    j4osc_init!(j4oscd::J4OsculatingPropagator, orb₀::KeplerianElements) -> Nothing

Initialize the J4 osculating orbit propagator structure `j4oscd` using the mean Keplerian
elements `orb₀`.

!!! warning

    The propagation constants `j4c::J4PropagatorConstants` in `j4oscd.j4d` will not be
    changed. Hence, they must be initialized.
"""
function j4osc_init!(j4oscd::J4OsculatingPropagator, orb₀::KeplerianElements)
    # Initialize the J4 propagator that will propagate the mean elements.
    j4_init!(j4oscd.j4d, orb₀)

    # Call the propagation one time to update the osculating elements.
    j4osc!(j4oscd, 0)

    return nothing
end

"""
    j4osc(Δt::Number, orb₀::KeplerianElements; kwargs...) -> SVector{3, T}, SVector{3, T}, J4OsculatingPropagator

Initialize the J4 osculating propagator structure using the input elements `orb₀` and
propagate the orbit until the time Δt [s].

!!! note

    The type used in the propagation will be the same as used to define the constants in the
    structure `j4c`.

# Keywords

- `j4c::J4PropagatorConstants{T}`: J4 orbit propagator constants (see
  [`J4PropagatorConstants`](@ref)). (**Default** = `j4c_egm2008`)

# Returns

- `SVector{3, T}`: Position vector [m] represented in the inertial frame at propagation
    instant.
- `SVector{3, T}`: Velocity vector [m / s] represented in the inertial frame at propagation
    instant.
- [`J4OsculatingPropagator`](@ref): Structure with the initialized parameters.

# Remarks

The inertial frame in which the output is represented depends on which frame it was used to
generate the orbit parameters. Notice that the perturbation theory requires an inertial
frame with true equator.
"""
function j4osc(Δt::Number, orb₀::KeplerianElements; j4c::J4PropagatorConstants = j4c_egm2008)
    j4oscd = j4osc_init(orb₀; j4c = j4c)
    r_i, v_i = j4osc!(j4oscd, Δt)
    return r_i, v_i, j4oscd
end

"""
    j4osc!(j4oscd::J4OsculatingPropagator{Tepoch, T}, t::Number) where {Tepoch, T} -> SVector{3, T}, SVector{3, T}

Propagate the orbit defined in `j4oscd` (see [`J4OsculatingPropagator`](@ref)) to `t` [s]
after the epoch of the input mean elements in `j4d`.

!!! note

    The internal values in `j4oscd` will be modified.

# Returns

- `SVector{3, T}`: Position vector [m] represented in the inertial frame at propagation
    instant.
- `SVector{3, T}`: Velocity vector [m / s] represented in the inertial frame at propagation
    instant.

# Remarks

The inertial frame in which the output is represented depends on which frame it was used to
generate the orbit parameters. Notice that the perturbation theory requires an inertial
frame with true equator.
"""
function j4osc!(j4oscd::J4OsculatingPropagator{Tepoch, T}, t::Number) where {Tepoch<:Number, T<:Number}
    # First, we need to propagate the mean elements since they are necessary to compute the
    # short-periodic perturbations.
    j4d = j4oscd.j4d
    j4!(j4d, t)

    # Unpack the propagator constants.
    j4c = j4d.j4c
    R₀  = j4c.R0
    μm  = j4c.μm
    J₂  = j4c.J2

    # Time from epoch to propagate the orbit.
    Δt = T(t)

    # Obtain the mean elements at this time instant.
    mean_orbk = j4d.orbk

    a_k  = mean_orbk.a
    e_k  = mean_orbk.e
    e_k² = e_k * e_k
    i_k  = mean_orbk.i
    Ω_k  = mean_orbk.Ω
    ω_k  = mean_orbk.ω
    f_k  = mean_orbk.f
    M_k  = true_to_mean_anomaly(e_k, f_k)
    p_k  = a_k * (1 - e_k²)
    p_k² = p_k * p_k
    u_k  = ω_k + f_k

    # Auxiliary variables to reduce the computational burden.
    KJ₂ = J₂ * R₀ * R₀

    sin_i_k, cos_i_k         = sincos(i_k)
    sin_f_k, cos_f_k         = sincos(f_k)
    sin_2u_k, cos_2u_k       = sincos(2u_k)
    sin_2ω_f_k, cos_2ω_f_k   = sincos(2ω_k + f_k)
    sin_2ω_3f_k, cos_2ω_3f_k = sincos(2ω_k + 3f_k)

    sin_i_k²  = sin_i_k * sin_i_k
    cos_i_k²  = cos_i_k * cos_i_k
    e_cos_f_k = e_k * cos_f_k
    e_sin_f_k = e_k * sin_f_k

    aux1 = 3cos_2u_k + 3e_k * cos_2ω_f_k + e_k * cos_2ω_3f_k
    aux2 = √(1 - e_k²)
    aux3 = 3cos_i_k² - 1
    aux4 = -aux3 / (1 + aux2)

    # Compute the short-periodic perturbations considering only the J2 gravitational term.
    δisp_k = +KJ₂ * sin_i_k * cos_i_k / (4p_k²) * aux1

    δpsp_k = +KJ₂ * sin_i_k² / (2p_k) * aux1

    δΩsp_k = -KJ₂ * cos_i_k / (4p_k²) * (
        6 * (f_k - M_k + e_sin_f_k) - 3sin_2u_k - 3e_k * sin_2ω_f_k - e_k * sin_2ω_3f_k
    )

    δrsp_k = -KJ₂ / (4p_k) * (
        aux3 * (2aux2 / (1 + e_cos_f_k) + e_cos_f_k / (1 + aux2) + 1) - sin_i_k² * cos_2u_k
    )

    δṙsp_k = +KJ₂ * √μm / (4 * √(p_k^5)) * (
        aux3 * e_sin_f_k * (aux2 + ((1 + e_cos_f_k)^2) / (1 + aux2)) -
        2sin_i_k² * (1 - e_cos_f_k)^2 * sin_2u_k
    )

    δusp_k = +KJ₂ / (8p_k²) * (
        (6 - 30cos_i_k²) * (f_k - M_k) +
        4e_sin_f_k * (1 - 6cos_i_k² + aux4) +
        aux4 * e_k² * sin(2f_k) +
        (5cos_i_k² - 2) * (2e_k) * sin_2ω_f_k +
        (7cos_i_k² - 1) * sin_2u_k +
        2cos_i_k² * e_k * sin_2ω_3f_k
    )

    r_k = p_k / (1 + e_cos_f_k)
    ṙ_k = √(μm / p_k) * e_sin_f_k

    r_osc_k = r_k + δrsp_k
    ṙ_osc_k = ṙ_k + δṙsp_k
    p_osc_k = p_k + δpsp_k

    A_k = p_osc_k / r_osc_k - 1
    B_k = √(p_osc_k / μm) * ṙ_osc_k

    e_osc_k² = A_k^2 + B_k^2
    e_osc_k  = √e_osc_k²
    a_osc_k  = p_osc_k / (1 - e_osc_k²)
    i_osc_k  = i_k + δisp_k
    Ω_osc_k  = Ω_k + δΩsp_k
    u_osc_k  = u_k + δusp_k
    f_osc_k  = atan(B_k, A_k)
    ω_osc_k  = u_osc_k - f_osc_k

    # Assemble the current osculating elements.
    orbk = KeplerianElements(
        j4d.orb₀.t + Tepoch(Δt) / 86400,
        a_osc_k,
        e_osc_k,
        i_osc_k,
        Ω_osc_k,
        ω_osc_k,
        f_osc_k
    )

    # Compute the position and velocity considering the osculating elements.
    r_i_k, v_i_k = kepler_to_rv(orbk)

    # Update the J4 orbit propagator structure.
    j4oscd.Δt   = T(t)
    j4oscd.orbk = orbk

    return r_i_k, v_i_k
end

"""
    fit_j4osc_mean_elements(vjd::AbstractVector{Tjd}, vr_i::AbstractVector{Tv}, vv_i::AbstractVector{Tv}; kwargs...) where {Tjd<:Number, Tv<:AbstractVector} -> KeplerianElements{Float64, Float64}, SMatrix{6, 6, Float64}

Fit a set of mean Keplerian elements for the J4 osculating orbit propagator using the
osculating elements represented by a set of position vectors `vr_i` [m] and a set of
velocity vectors `vv_i` [m / s] represented in an inertial reference frame at instants in
the array `vjd` [Julian Day].

!!! note

    This algorithm version will allocate a new J4 osculating propagator with the default
    constants `j4c_egm2008`. If another set of constants are required, use the function
    [`fit_j4osc_mean_elements!`](@ref) instead.

# Keywords

- `atol::Number`: Tolerance for the residue absolute value. If the residue is lower than
    `atol` at any iteration, the computation loop stops.
    (**Default** = 2e-4)
- `rtol::Number`: Tolerance for the relative difference between the residues. If the
    relative difference between the residues in two consecutive iterations is lower than
    `rtol`, the computation loop stops.
    (**Default** = 2e-4)
- `initial_guess::Union{Nothing, KeplerianElements}`: Initial guess for the mean elements
    fitting process. If it is `nothing`, the algorithm will obtain an initial estimate from
    the osculating elements in `vr_i` and `vv_i`.
    (**Default** = nothing)
- `jacobian_perturbation::Number`: Initial state perturbation to compute the
    finite-difference when calculating the Jacobian matrix.
    (**Default** = 1e-3)
- `jacobian_perturbation_tol::Number`: Tolerance to accept the perturbation when calculating
    the Jacobian matrix. If the computed perturbation is lower than
    `jacobian_perturbation_tol`, we increase it until it absolute value is higher than
    `jacobian_perturbation_tol`.
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

julia> orb, P = fit_j4osc_mean_elements(vjd, vr_i, vv_i)
ACTION:   Fitting the mean elements for the J4 osculating propagator.
           Iteration        Position RMSE        Velocity RMSE           Total RMSE       RMSE Variation
                                     [km]             [km / s]                  [ ]
PROGRESS:          4          1.68852e-05           0.00259716              2.59722          6.52143e-10 %

(KeplerianElements{Float64, Float64}: Epoch = 2.46003e6 (2023-03-24T16:33:40.388), [0.999942775712293 -0.004901178682591695 … -0.00012021347649120565 -0.0001155014732440418; -0.004901178674245216 1.0007325984965127 … 0.0032661543862499806 3.8581997167487534e-5; … ; -0.00012021347646254671 0.003266154386250846 … 2.204822797002821e-5 5.3982040158043245e-8; -0.00011550147323899825 3.8581997169916946e-5 … 5.398204016417114e-8 2.1657607378172235e-5])

julia> orb
KeplerianElements{Float64, Float64}:
           Epoch :    2.46003e6 (2023-03-24T16:33:40.388)
 Semi-major axis : 7135.79       km
    Eccentricity :    0.00135314
     Inclination :   98.4304     °
            RAAN :  162.113      °
 Arg. of Perigee :   64.9687     °
    True Anomaly :  313.042      °
```
"""
function fit_j4osc_mean_elements(
    vjd::AbstractVector{Tjd},
    vr_i::AbstractVector{Tv},
    vv_i::AbstractVector{Tv};
    kwargs...
) where {Tjd<:Number, Tv<:AbstractVector}
    # Allocate the J4 propagator structure that will propagate the mean elements.
    j4d = J4Propagator{Float64, Float64}()

    # Assign the constants, which are used in the initialization.
    j4d.j4c = j4c_egm2008

    # Allocate the J4 osculating propagator structure.
    j4oscd = J4OsculatingPropagator{Float64, Float64}()
    j4oscd.j4d = j4d

    return fit_j4osc_mean_elements!(j4oscd, vjd, vr_i, vv_i; kwargs...)
end

"""
    fit_j4osc_mean_elements!(j4oscd::J4OsculatingPropagator{Tepoch, T}, vjd::AbstractVector{Tjd}, vr_i::AbstractVector{Tv}, vv_i::AbstractVector{Tv}; kwargs...) where {T<:Number, Tepoch<:Number, Tjd<:Number, Tv<:AbstractVector} -> KeplerianElements{Tepoch, T}, SMatrix{6, 6, T}

Fit a set of mean Keplerian elements for the J4 osculating orbit propagator `j4oscd` using
the osculating elements represented by a set of position vectors `vr_i` [m] and a set of
velocity vectors `vv_i` [m / s] represented in an inertial reference frame at instants in
the array `vjd` [Julian Day].

!!! note

    The J4 osculating orbit propagator `j4oscd` will be initialized with the Keplerian
    elements returned by the function.

# Keywords

- `atol::Number`: Tolerance for the residue absolute value. If the residue is lower than
    `atol` at any iteration, the computation loop stops.
    (**Default** = 2e-4)
- `rtol::Number`: Tolerance for the relative difference between the residues. If the
    relative difference between the residues in two consecutive iterations is lower than
    `rtol`, the computation loop stops.
    (**Default** = 2e-4)
- `initial_guess::Union{Nothing, KeplerianElements}`: Initial guess for the mean elements
    fitting process. If it is `nothing`, the algorithm will obtain an initial estimate from
    the osculating elements in `vr_i` and `vv_i`.
    (**Default** = nothing)
- `jacobian_perturbation::Number`: Initial state perturbation to compute the
    finite-difference when calculating the Jacobian matrix.
    (**Default** = 1e-3)
- `jacobian_perturbation_tol::Number`: Tolerance to accept the perturbation when calculating
    the Jacobian matrix. If the computed perturbation is lower than
    `jacobian_perturbation_tol`, we increase it until it absolute value is higher than
    `jacobian_perturbation_tol`.
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
# Allocate a new J4 osculating orbit propagator using a dummy Keplerian elements.
julia> j4oscd = j4osc_init(KeplerianElements{Float64, Float64}(0, 7000e3, 0, 0, 0, 0, 0));

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

julia> orb, P = fit_j4osc_mean_elements!(j4oscd, vjd, vr_i, vv_i)
ACTION:   Fitting the mean elements for the J4 osculating propagator.
           Iteration        Position RMSE        Velocity RMSE           Total RMSE       RMSE Variation
                                     [km]             [km / s]                  [ ]
PROGRESS:          4          1.68852e-05           0.00259716              2.59722          6.52143e-10 %

(KeplerianElements{Float64, Float64}: Epoch = 2.46003e6 (2023-03-24T16:33:40.388), [0.999942775712293 -0.004901178682591695 … -0.00012021347649120565 -0.0001155014732440418; -0.004901178674245216 1.0007325984965127 … 0.0032661543862499806 3.8581997167487534e-5; … ; -0.00012021347646254671 0.003266154386250846 … 2.204822797002821e-5 5.3982040158043245e-8; -0.00011550147323899825 3.8581997169916946e-5 … 5.398204016417114e-8 2.1657607378172235e-5])

julia> orb
KeplerianElements{Float64, Float64}:
           Epoch :    2.46003e6 (2023-03-24T16:33:40.388)
 Semi-major axis : 7135.79       km
    Eccentricity :    0.00135314
     Inclination :   98.4304     °
            RAAN :  162.113      °
 Arg. of Perigee :   64.9687     °
    True Anomaly :  313.042      °
```
"""
function fit_j4osc_mean_elements!(
    j4oscd::J4OsculatingPropagator{Tepoch, T},
    vjd::AbstractVector{Tjd},
    vr_i::AbstractVector{Tv},
    vv_i::AbstractVector{Tv};
    atol::Number                                     = 2e-4,
    rtol::Number                                     = 2e-4,
    initial_guess::Union{Nothing, KeplerianElements} = nothing,
    jacobian_perturbation::Number                    = 1e-3,
    jacobian_perturbation_tol::Number                = 1e-7,
    max_iterations::Int                              = 50,
    mean_elements_epoch::Number                      = vjd[end],
    verbose::Bool                                    = true,
    weight_vector::AbstractVector                    = @SVector(ones(Bool, 6)),
) where {T<:Number, Tepoch<:Number, Tjd<:Number, Tv<:AbstractVector}
    # Number of available measurements.
    num_measurements = length(vjd)

    # Check the inputs.
    length(vr_i) != num_measurements &&
        throw(ArgumentError("The number of elements in `vjd` and `vr_i` must be the same."))

    length(vv_i) != num_measurements &&
        throw(ArgumentError("The number of elements in `vjd` and `vv_i` must be the same."))

    if length(weight_vector) != 6
        throw(ArgumentError("The weight vector must have 6 elements."))
    end

    # Check if stdout supports colors.
    has_color = get(stdout, :color, false)::Bool
    cd = has_color ? _D : ""
    cb = has_color ? _B : ""
    cy = has_color ? _Y : ""

    # Assemble the weight matrix.
    W = Diagonal(
        @SVector T[
            weight_vector[1],
            weight_vector[2],
            weight_vector[3],
            weight_vector[4],
            weight_vector[5],
            weight_vector[6]
        ]
    )

    # Initial guess of the mean elements.
    #
    # NOTE: x₁ is the previous estimate and x₂ is the current estimate.
    if !isnothing(initial_guess)
        epoch = mean_elements_epoch

        # First, we need to update the mean elements to the desired epoch.
        verbose && println("$(cy)ACTION:$(cd)   Updating the epoch of the initial mean elements guess to match the desired one.")
        orb = update_j4osc_mean_elements_epoch!(j4oscd, initial_guess, epoch)

        r_i, v_i = kepler_to_rv(orb)
        x₁ = SVector{6, T}(r_i..., v_i...)
    else
        # In this case, we must find the closest osculating vector to the desired epoch.
        id = firstindex(vjd)
        v  = abs(vjd[1] - mean_elements_epoch)

        for k in eachindex(vjd)
            vk = abs(vjd[k] - mean_elements_epoch)
            if vk < v
                id = k
                v  = vk
            end
        end

        epoch = T(vjd[id])
        r_i   = vr_i[id]
        v_i   = vv_i[id]
        x₁    = SVector{6, T}(r_i[1], r_i[2], r_i[3], v_i[1], v_i[2], v_i[3],)
    end

    x₂ = x₁

    # Number of states in the input vector.
    num_states = 6

    # Covariance matrix.
    P = SMatrix{num_states, num_states, T}(I)

    # Variable to store the last residue.
    local σ_i_₁

    # Variable to store how many iterations the residue increased. This is used to account
    # for divergence.
    Δd = 0

    # Header.
    if verbose
        println("$(cy)ACTION:$(cd)   Fitting the mean elements for the J4 osculating propagator.")
        @printf("          %s%10s %20s %20s %20s %20s%s\n", cy, "Iteration", "Position RMSE", "Velocity RMSE", "Total RMSE", "RMSE Variation", cd)
        @printf("          %s%10s %20s %20s %20s %20s%s\n", cb, "", "[km]", "[km / s]", "[ ]", "", cd)
        println()
    end

    # We need a reference to the covariance inverse because we will invert it and return
    # after the iterations.
    local ΣJ′WJ

    # Loop until the maximum allowed iteration.
    @inbounds @views for it in 1:max_iterations
        x₁ = x₂

        # Variables to store the summations to compute the least square fitting algorithm.
        ΣJ′WJ = @SMatrix zeros(T, num_states, num_states)
        ΣJ′Wb = @SVector zeros(T, num_states)

        # Variable to store the RMS errors in this iteration.
        σ_i  = T(0)
        σp_i = T(0)
        σv_i = T(0)

        for k in 1:num_measurements
            # Obtain the measured ephemerides.
            y = vcat(vr_i[k - 1 + begin], vv_i[k - 1 + begin])

            # Initialize the SGP4 with the current estimated mean elements.
            orb = rv_to_kepler(x₁[1:3], x₁[4:6], epoch)
            j4osc_init!(j4oscd, orb)

            # Obtain the propagation time for this measurement.
            Δt = (vjd[k - 1 + begin] - epoch) * 86400

            # Propagate the orbit.
            r̂_i, v̂_i = j4osc!(j4oscd, Δt)
            ŷ = vcat(r̂_i, v̂_i)

            # Compute the residue.
            b = y - ŷ

            # Compute the Jacobian in-place.
            J = _j4osc_jacobian(
                j4oscd,
                Δt,
                x₁,
                ŷ;
                perturbation     = jacobian_perturbation,
                perturbation_tol = jacobian_perturbation_tol
            )

            # Accumulation.
            ΣJ′WJ += J' * W * J
            ΣJ′Wb += J' * W * b
            σ_i   += b' * W * b
            σp_i  += dot(b[1:3], b[1:3])
            σv_i  += dot(b[4:6], b[4:6])
        end

        # Normalize and compute the RMS errors.
        σ_i  = √(σ_i  / num_measurements)
        σp_i = √(σp_i / num_measurements)
        σv_i = √(σv_i / num_measurements)

        # Update the estimate.
        δx = ΣJ′WJ \ ΣJ′Wb

        # Limit the correction to avoid divergence.
        for i in 1:num_states
            threshold = T(0.1)
            if abs(δx[i] / x₁[i]) > threshold
                δx = setindex(δx, threshold * abs(x₁[i]) * sign(δx[i]), i)
            end
        end

        x₂ = x₁ + δx

        # We cannot compute the RMSE variation in the first iteration.
        if it == 1
            verbose &&
                @printf("\x1b[A\x1b[2K\r%sPROGRESS:%s %10d %20g %20g %20g %20s\n", cb, cd, it, σp_i / 1000, σv_i / 1000, σ_i, "---")

        else
            # Compute the RMSE variation.
            Δσ = (σ_i - σ_i_₁) / σ_i_₁

            verbose &&
                @printf("\x1b[A\x1b[2K\r%sPROGRESS:%s %10d %20g %20g %20g %20g %%\n", cb, cd, it, σp_i / 1000, σv_i / 1000, σ_i, 100 * Δσ)

            # Check if the RMSE is increasing.
            if σ_i < σ_i_₁
                Δd = 0
            else
                Δd += 1
            end

            # If the RMSE increased by three iterations and its value is higher than 5e11,
            # we abort because the iterations are diverging.
            ((Δd ≥ 3) && (σ_i > 5e11)) && error("The iterations diverged!")

            # Check if the condition to stop has been reached.
            ((abs(Δσ) < rtol) || (σ_i < atol) || (it ≥ max_iterations)) && break
        end

        σ_i_₁ = σ_i
    end

    verbose && println()

    # Obtain the mean elements.
    orb = @views rv_to_kepler(x₂[1:3], x₂[4:6], epoch)

    # Update the epoch of the fitted mean elements to match the desired one.
    if abs(epoch - mean_elements_epoch) > 0.001 / 86400
        verbose &&
            println("$(cy)ACTION:$(cd)   Updating the epoch of the fitted mean elements to match the desired one.")
        orb = update_j4osc_mean_elements_epoch!(j4oscd, orb, mean_elements_epoch)
    end

    # Initialize the propagator with the mean elements.
    j4osc_init!(j4oscd, orb)

    # Compute the final covariance.
    P = pinv(ΣJ′WJ)

    # Return the mean elements and the covariance.
    return orb, P
end

"""
    update_j4osc_mean_elements_epoch(orb::KeplerianElements, new_epoch::Union{Number, DateTime}) -> KepleriranElements

Update the epoch of the mean elements `orb` using a J4 osculating orbit propagator to
`new_epoch`, which can be represented by a Julian Day or a `DateTime`.

!!! note

    This algorithm version will allocate a new J4 osculating propagator with the default
    constants `j4c_egm2008`. If another set of constants are required, use the function
    [`update_j4osc_mean_elements_epoch!`](@ref) instead.

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

julia> update_j4osc_mean_elements_epoch(orb, DateTime("2023-01-02"))
KeplerianElements{Float64, Float64}:
           Epoch :    2.45995e6 (2023-01-02T00:00:00)
 Semi-major axis : 7190.98     km
    Eccentricity :    0.001111
     Inclination :   98.405    °
            RAAN :   90.9555   °
 Arg. of Perigee :  197.079    °
    True Anomaly :  127.293    °
```
"""
function update_j4osc_mean_elements_epoch(
    orb::KeplerianElements{Tepoch, T},
    new_epoch::Union{Number, DateTime}
) where {T<:Number, Tepoch<:Number}
    # Allocate the J4 propagator structure that will propagate the mean elements.
    j4d = J4Propagator{Tepoch, T}()

    # Assign the constants, which are used in the initialization.
    j4d.j4c = j4c_egm2008

    # Allocate the J4 osculating propagator structure.
    j4oscd = J4OsculatingPropagator{Tepoch, T}()
    j4oscd.j4d = j4d

    return update_j4osc_mean_elements_epoch!(j4oscd, orb, new_epoch)
end

"""
    update_j4osc_mean_elements_epoch!(j4oscd::J4OsculatingPropagator, orb::KeplerianElements, new_epoch::Union{Number, DateTime}) -> KepleriranElements

Update the epoch of the mean elements `orb` using the propagator `j4oscd` to `new_epoch`,
which can be represented by a Julian Day or a `DateTime`.

!!! note

    The J4 osculating orbit propagator `j4oscd` will be initialized with the Keplerian
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

# Allocate a new J4 osculating orbit propagator using the created Keplerian elements. Notice
# that any set of Keplerian elements can be used here.
julia> j4oscd = j4osc_init(orb);

julia> update_j4osc_mean_elements_epoch!(j4oscd, orb, DateTime("2023-01-02"))
KeplerianElements{Float64, Float64}:
           Epoch :    2.45995e6 (2023-01-02T00:00:00)
 Semi-major axis : 7190.98     km
    Eccentricity :    0.001111
     Inclination :   98.405    °
            RAAN :   90.9555   °
 Arg. of Perigee :  197.079    °
    True Anomaly :  127.293    °
```
"""
function update_j4osc_mean_elements_epoch!(
    j4oscd::J4OsculatingPropagator,
    orb::KeplerianElements,
    new_epoch::DateTime
)
    dt = datetime2julian(new_epoch)
    return update_j4osc_mean_elements_epoch!(j4oscd, orb, dt)
end

function update_j4osc_mean_elements_epoch!(
    j4oscd::J4OsculatingPropagator,
    orb::KeplerianElements,
    new_epoch::Number
)
    # First, we need to initialize the J4 osculating propagator with the mean elements.
    j4osc_init!(j4oscd, orb)

    # Now, we just need to propagate the orbit to the desired instant and obtain the mean
    # elements from the J4 propagator structure inside.
    Δt = (new_epoch - j4oscd.j4d.orb₀.t) * 86400
    j4osc!(j4oscd, Δt)
    orb = j4oscd.j4d.orbk

    # Finally, we initialize the propagator with the new set of mean elements.
    j4osc_init!(j4oscd, orb)

    return orb
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

"""
    _j4osc_jacobian(j4oscd::J4OsculatingPropagator{Tepoch, T}, Δt::Number, x₁::SVector{6, T}, y₁::SVector{6, T}; kwargs...)) where {T<:Number, Tepoch<:Number} -> SMatrix{6, 6, T}

Compute the J4 osculating orbit propagator Jacobian by finite-differences using the
propagator `j4oscd` at instant `Δt` considering the input mean elements `x₁` that must
provide the output vector `y₁`. Hence:

        ∂j4osc(x, Δt) │
    J = ───────────── │
              ∂x      │ x = x₁

# Keywords

- `perturbation::T`: Initial state perturbation to compute the finite-difference:
    `Δx = x * perturbation`.
    (**Default** = 1e-3)
- `perturbation_tol::T`: Tolerance to accept the perturbation. If the computed perturbation
    is lower than `perturbation_tol`, we increase it until it absolute value is higher than
    `perturbation_tol`.
    (**Default** = 1e-7)
"""
function _j4osc_jacobian(
    j4oscd::J4OsculatingPropagator{Tepoch, T},
    Δt::Number,
    x₁::SVector{6, T},
    y₁::SVector{6, T};
    perturbation::Number = T(1e-3),
    perturbation_tol::Number = T(1e-7)
) where {T<:Number, Tepoch<:Number}

    # Allocate the `MMatrix` that will have the Jacobian.
    J = MMatrix{6, 6, T}(undef)

    # Auxiliary variables.
    x₂ = x₁

    @inbounds @views for j in 1:6
        # State that will be perturbed.
        α = x₂[j]

        # Obtain the perturbation, taking care to avoid small values.
        ϵ = α * T(perturbation)

        for _ in 1:5
            abs(ϵ) > perturbation_tol && break
            ϵ *= T(1.4)
        end

        # Avoid division by zero in cases that α is very small. In this situation, we force
        # `|ϵ| = perturbation_tol`.
        if abs(ϵ) < perturbation_tol
            ϵ = signbit(α) ? -perturbation_tol : perturbation_tol
        end

        α += ϵ

        # Modify the perturbed state.
        x₂ = setindex(x₂, α, j)

        # Obtain the Jacobian by finite differentiation.
        orb = rv_to_kepler(x₂[1:3], x₂[4:6], j4oscd.j4d.orb₀.t)
        j4osc_init!(j4oscd, orb)
        r_i, v_i = j4osc!(j4oscd, Δt)
        y₂ = @SVector [
            r_i[1],
            r_i[2],
            r_i[3],
            v_i[1],
            v_i[2],
            v_i[3],
        ]

        J[:, j] .= (y₂ .- y₁) ./ ϵ

        # Restore the value of the perturbed state for the next iteration.
        x₂ = setindex(x₂, x₁[j], j)
    end

    # Convert the Jacobian to a static matrix to avoid allocations.
    Js = SMatrix{6, 6, T}(J)

    return Js
end
