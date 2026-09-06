## Description #############################################################################
#
# Functions to print the structures of the orbit propagators implemented in this package.
#
# The layout follows the orbit representations of SatelliteToolboxBase.jl: the compact form
# prints the type with its parameters and the epoch, whereas the rich form prints one field
# per line, aligned at the decimal point. The `OrbitPropagator` wrappers of the API delegate
# to those representations in `src/api/Propagators.jl`.
#
############################################################################################

############################################################################################
#                                      Compact Format                                      #
############################################################################################

function Base.show(io::IO, pd::PropagatorData)
    name = Propagators._type_name(pd)

    if !_is_initialized(pd)
        print(io, name, " (not initialized)")
        return nothing
    end

    SatelliteToolboxBase.print_compact(io, name, _initial_epoch(pd))

    return nothing
end

############################################################################################
#                                       Rich Format                                        #
############################################################################################

function Base.show(io::IO, ::MIME"text/plain", pd::PropagatorData)
    name = Propagators._type_name(pd)

    if !_is_initialized(pd)
        println(io, name, ":")
        SatelliteToolboxBase.print_field(io, " Status : ", "not initialized")
        return nothing
    end

    labels, values, units = _show_rows(io, pd)

    SatelliteToolboxBase.print_elements(
        io,
        name,
        _initial_epoch(pd),
        (labels..., "Last propagation"),
        (values..., SatelliteToolboxBase.compact_string(io, pd.Δt)),
        (units..., "s"),
    )

    return nothing
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

"""
    _initial_epoch(pd::PropagatorData) -> Number

Return the epoch [Julian Day] of the initial mean elements stored in the initialized
propagator structure `pd`.
"""
_initial_epoch(pd::J2Propagator)           = pd.orb₀.epoch
_initial_epoch(pd::J2OsculatingPropagator) = _initial_epoch(pd.j2d)
_initial_epoch(pd::J4Propagator)           = pd.orb₀.epoch
_initial_epoch(pd::J4OsculatingPropagator) = _initial_epoch(pd.j4d)
_initial_epoch(pd::TwoBodyPropagator)      = pd.orb₀.epoch

"""
    _is_initialized(pd::PropagatorData) -> Bool

Return whether the propagator structure `pd` has been initialized. The osculating
propagators created with the empty constructor have an undefined field `j2d` or `j4d` until
they are initialized. The other structures only have `isbits` fields, which always hold a
value, so they cannot be distinguished from initialized ones and this function returns
`true`.
"""
_is_initialized(pd::J2Propagator)           = true
_is_initialized(pd::J2OsculatingPropagator) = isdefined(pd, :j2d)
_is_initialized(pd::J4Propagator)           = true
_is_initialized(pd::J4OsculatingPropagator) = isdefined(pd, :j4d)
_is_initialized(pd::TwoBodyPropagator)      = true

"""
    _mean_elements_rows(orb₀::KeplerianElements{MeanAnomaly}) -> Tuple, Tuple, Tuple

Return the labels, the values, and the units of the rows that print the initial mean
elements `orb₀` in the rich representation of a propagator structure. The semi-major axis
is printed in kilometers and the angles in degrees, all with 8 decimal digits.
"""
function _mean_elements_rows(orb₀::KeplerianElements{MeanAnomaly})
    labels = (
        "Semi-major axis",
        "Eccentricity",
        "Inclination",
        "RAAN",
        "Arg. of perigee",
        "Mean anomaly",
    )

    values = (
        _show_number(orb₀.semi_major_axis / 1000),
        _show_number(orb₀.eccentricity),
        _show_number(rad2deg(orb₀.inclination)),
        _show_number(rad2deg(orb₀.raan)),
        _show_number(rad2deg(orb₀.argument_of_periapsis)),
        _show_number(rad2deg(orb₀.anomaly)),
    )

    units = ("km", "", "°", "°", "°", "°")

    return labels, values, units
end

"""
    _secular_rates_rows(n̄::Number, ∂Ω::Number, ∂ω::Number) -> Tuple, Tuple, Tuple

Return the labels, the values, and the units of the rows that print the perturbed mean
motion `n̄` [rad / s], the RAAN rate `∂Ω` [rad / s], and the argument of perigee rate `∂ω`
[rad / s] in the rich representation of a propagator structure. The mean motion is printed
in revolutions per day and the rates in degrees per day.
"""
function _secular_rates_rows(n̄::Number, ∂Ω::Number, ∂ω::Number)
    labels = ("Mean motion", "RAAN rate", "Arg. of perigee rate")

    values = (
        _show_number(86400 * n̄ / 2π),
        _show_number(86400 * rad2deg(∂Ω)),
        _show_number(86400 * rad2deg(∂ω)),
    )

    units = ("rev / day", "° / day", "° / day")

    return labels, values, units
end

"""
    _show_number(x::Number) -> String

Format the number `x` to be printed with 8 decimal digits if it is a floating-point number.
Otherwise, e.g. for dual numbers, it is printed as is.
"""
_show_number(x::AbstractFloat) = @sprintf("%.8f", x)
_show_number(x::Number) = string(x)

"""
    _show_rows(io::IO, pd::PropagatorData) -> Tuple, Tuple, Tuple

Return the labels, the values, and the units of the rows printed by the rich representation
of the initialized propagator structure `pd`, excluding the epoch and the last propagation
instant, which are printed by the caller. The constants are printed using the `:compact`
property of `io`.
"""
function _show_rows(io::IO, pd::J2Propagator)
    j2c = pd.j2c

    labels = ("R₀", "μm", "J₂")

    values = (
        SatelliteToolboxBase.compact_string(io, j2c.R0 / 1000),
        SatelliteToolboxBase.compact_string(io, j2c.μm),
        SatelliteToolboxBase.compact_string(io, j2c.J2),
    )

    units = ("km", "rad / s", "")

    orb_labels, orb_values, orb_units = _mean_elements_rows(pd.orb₀)
    rate_labels, rate_values, rate_units = _secular_rates_rows(pd.n̄, pd.∂Ω, pd.∂ω)

    return (
        (labels..., orb_labels..., rate_labels...),
        (values..., orb_values..., rate_values...),
        (units..., orb_units..., rate_units...),
    )
end

function _show_rows(io::IO, pd::J4Propagator)
    j4c = pd.j4c

    labels = ("R₀", "μm", "J₂", "J₄")

    values = (
        SatelliteToolboxBase.compact_string(io, j4c.R0 / 1000),
        SatelliteToolboxBase.compact_string(io, j4c.μm),
        SatelliteToolboxBase.compact_string(io, j4c.J2),
        SatelliteToolboxBase.compact_string(io, j4c.J4),
    )

    units = ("km", "rad / s", "", "")

    orb_labels, orb_values, orb_units = _mean_elements_rows(pd.orb₀)
    rate_labels, rate_values, rate_units = _secular_rates_rows(pd.n̄, pd.∂Ω, pd.∂ω)

    return (
        (labels..., orb_labels..., rate_labels...),
        (values..., orb_values..., rate_values...),
        (units..., orb_units..., rate_units...),
    )
end

function _show_rows(io::IO, pd::TwoBodyPropagator)
    labels = ("μ",)
    values = (SatelliteToolboxBase.compact_string(io, pd.μ),)
    units  = ("m³ / s²",)

    orb_labels, orb_values, orb_units = _mean_elements_rows(pd.orb₀)

    return (
        (labels..., orb_labels..., "Mean motion"),
        (values..., orb_values..., _show_number(86400 * pd.n₀ / 2π)),
        (units..., orb_units..., "rev / day"),
    )
end

# The osculating propagators print the same rows as the propagators of the mean elements
# they wrap, since the short-period corrections do not add any parameter.
_show_rows(io::IO, pd::J2OsculatingPropagator) = _show_rows(io, pd.j2d)
_show_rows(io::IO, pd::J4OsculatingPropagator) = _show_rows(io, pd.j4d)
