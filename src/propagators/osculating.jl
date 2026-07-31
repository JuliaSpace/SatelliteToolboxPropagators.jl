## Description #############################################################################
#
# Short-period corrections introduced by the J2 gravitational term.
#
# Both the J2 and the J4 osculating propagators apply the same corrections to the mean
# elements, since the J4 propagator only adds secular terms. Hence, this algorithm is
# implemented once here and used by both.
#
## References ##############################################################################
#
# [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. 4th ed.
#     Microcosm Press, Hawthorn, CA, USA.
#
############################################################################################

# Compute the osculating elements at the epoch `epoch` given the mean elements `mean_orbk`,
# the mean anomaly `M_k` [rad] associated with them, and the propagator constants `R₀` [m],
# `μm` [rad / s], and `J₂`. The correction considers only the J2 gravitational term, as
# described in [1, p. 708-710].
function _osculating_elements(
    mean_orbk::KeplerianElements{Tepoch, T},
    M_k::Number,
    epoch::Number,
    R₀::Number,
    μm::Number,
    J₂::Number,
) where {Tepoch <: Number, T <: Number}
    a_k  = mean_orbk.a
    e_k  = mean_orbk.e
    e_k² = e_k * e_k
    i_k  = mean_orbk.i
    Ω_k  = mean_orbk.Ω
    ω_k  = mean_orbk.ω
    f_k  = mean_orbk.f
    p_k  = a_k * (1 - e_k²)
    p_k² = p_k * p_k
    u_k  = ω_k + f_k

    # Auxiliary variable to reduce the computational burden.
    KJ₂ = J₂ * R₀ * R₀ / 4

    sin_i_k, cos_i_k = sincos(i_k)
    sin_f_k, cos_f_k = sincos(f_k)
    sin_2u_k, cos_2u_k = sincos(2u_k)

    # Since `u = ω + f`, we have `2ω + f == 2u - f` and `2ω + 3f == 2u + f`. Hence, we can
    # obtain the following sines and cosines from the ones already computed.
    sin_2ω_f_k  = sin_2u_k * cos_f_k - cos_2u_k * sin_f_k
    cos_2ω_f_k  = cos_2u_k * cos_f_k + sin_2u_k * sin_f_k
    sin_2ω_3f_k = sin_2u_k * cos_f_k + cos_2u_k * sin_f_k
    cos_2ω_3f_k = cos_2u_k * cos_f_k - sin_2u_k * sin_f_k

    sin_i_k²  = sin_i_k * sin_i_k
    cos_i_k²  = cos_i_k * cos_i_k
    e_cos_f_k = e_k * cos_f_k
    e_sin_f_k = e_k * sin_f_k

    # Both `f_k` and `M_k` are wrapped to [0, 2π). Hence, a rounding difference near the
    # boundary can make this difference jump by 2π, which would otherwise inject a
    # discontinuity of about 0.5° in the corrections below.
    Δf_M = rem2pi(f_k - M_k, RoundNearest)

    # Auxiliary variables to reduce the computational burden.
    k₁  = 3cos_2u_k + 3e_k * cos_2ω_f_k + e_k * cos_2ω_3f_k
    k₂  = √(1 - e_k²)
    k₃  = 3cos_i_k² - 1
    k₄  = -k₃ / (1 + k₂)
    k₅  = 1 + e_cos_f_k
    k₅² = k₅ * k₅

    # Notice that `μm` is √(GM / R0³) and not the standard gravitational parameter. Hence,
    # `ṙ_k` and `ṙ_osc_k` are not the physical radial rates. This is not a problem because
    # `μm` cancels out exactly when computing `B_k`, which is what `kepler_to_rv` requires.
    # However, all the three occurrences must be kept consistent.
    sqrt_μm_p = √(μm / p_k)

    # Compute the short-periodic perturbations considering only the J2 gravitational term.
    δisp_k = +KJ₂ * sin_i_k * cos_i_k / p_k² * k₁

    δpsp_k = +2KJ₂ * sin_i_k² / p_k * k₁

    δΩsp_k =
        -KJ₂ * cos_i_k / p_k² *
        (6 * (Δf_M + e_sin_f_k) - 3sin_2u_k - 3e_k * sin_2ω_f_k - e_k * sin_2ω_3f_k)

    δrsp_k = -KJ₂ / p_k * (k₃ * (2k₂ / k₅ + e_cos_f_k / (1 + k₂) + 1) - sin_i_k² * cos_2u_k)

    δṙsp_k = +KJ₂ * sqrt_μm_p / p_k² * (
        k₃ * e_sin_f_k * (k₂ + k₅² / (1 + k₂)) - 2sin_i_k² * k₅² * sin_2u_k
    )

    δusp_k =
        +KJ₂ / (2p_k²) * (
            (6 - 30cos_i_k²) * Δf_M +
            4e_sin_f_k * (1 - 6cos_i_k² + k₄) +
            k₄ * e_k² * 2sin_f_k * cos_f_k +
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
    Ω_osc_k  = mod(Ω_k + δΩsp_k, T(2π))
    u_osc_k  = u_k + δusp_k
    f_osc_k  = mod(atan(B_k, A_k), T(2π))
    ω_osc_k  = mod(u_osc_k - f_osc_k, T(2π))

    # Assemble the current osculating elements.
    return KeplerianElements(epoch, a_osc_k, e_osc_k, i_osc_k, Ω_osc_k, ω_osc_k, f_osc_k)
end
