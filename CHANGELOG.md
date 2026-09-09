SatelliteToolboxPropagator.jl Changelog
=======================================

Version 2.0.0
-------------

- ![BREAKING][badge-breaking] Require **SatelliteToolboxBase.jl** v2.1. The propagators store
  their mean elements as `KeplerianElements{MeanAnomaly}`, and the osculating propagators
  store the osculating elements as `KeplerianElements{TrueAnomaly}`. The fields `orb₀` and
  `orbk` of the propagator structures changed their types accordingly.
- ![BREAKING][badge-breaking] Rename the exported constants to `SCREAMING_SNAKE_CASE`:
  `J2C_EGM2008`, `J2C_EGM1996`, `J2C_JGM02`, `J2C_JGM03`, `J4C_EGM2008`, `J4C_EGM1996`,
  `J4C_JGM02`, `J4C_JGM03`, their `_F32` variants, `TBC_M0`, and `TBC_M0_F32`. The lowercase
  names were removed.
- ![BREAKING][badge-breaking] Every mean elements output now stores the mean anomaly. The
  functions that fit mean elements and the functions that update their epoch return
  `KeplerianElements{MeanAnomaly}`, and `Propagators.mean_elements` returns
  `KeplerianElements{MeanAnomaly}` for every propagator, including SGP4, which returned the
  true anomaly.
- ![BREAKING][badge-breaking] Require **SatelliteToolboxSgp4.jl** v3. The SGP4 constants
  are now `SGP4C_WGS84` and `SGP4C_WGS72`, and their `Float32` variants are obtained with
  `Sgp4Constants{Float32}(SGP4C_WGS84)`. An uninitialized SGP4 propagator is created with
  `Sgp4Propagator(SGP4C_WGS84)`, which also sets the deep space data.
- ![BREAKING][badge-breaking] The SGP4 fitting functions `Propagators.fit_mean_elements`
  and `Propagators.fit_mean_elements!` return an `OrbitMeanElementsMessage` by default. The
  representation is selected by the new sink argument, which can be `TLE` or
  `OrbitMeanElementsMessage`, placed before the propagator tag in the non-mutating function
  and after the measurements in the mutating one. The metadata of the output is provided by
  the keyword `template` instead of the six TLE-specific keywords, and a diverging fit
  throws `Sgp4FitDivergenceError` instead of an `ErrorException`.
- ![BREAKING][badge-breaking] Every function that fits the mean elements
  (`fit_*_mean_elements`, `fit_*_mean_elements!`, `Propagators.fit_mean_elements`, and
  `Propagators.fit_mean_elements!`) returns a third value: a `NamedTuple` with the
  statistics of the least-square algorithm, whose fields are `converged`, `iterations`,
  `position_rmse` [m], `velocity_rmse` [m / s], and `total_rmse`, as in
  **SatelliteToolboxSgp4.jl** v3. The SGP4 statistics are converted to SI units.
- ![BREAKING][badge-breaking] The fitting functions of the analytical propagators throw the
  new exported exception `MeanElementsFitDivergenceError`, which stores the iteration and
  the residue, when the least-square iterations diverge, instead of an `ErrorException`,
  and validate `max_iterations` with an `ArgumentError`, as **SatelliteToolboxSgp4.jl** v3
  does.
- ![Feature][badge-feature] Initialize the SGP4 orbit propagator with an Orbit
  Mean-Elements Message (OMM) using `Propagators.init(Val(:SGP4), omm)` and
  `Propagators.init!(orbp, omm)`, where `omm` is an `OrbitMeanElementsMessage` from
  **SatelliteToolboxOrbitDataMessages.jl**, which is re-exported through
  **SatelliteToolboxSgp4.jl**.
- ![Feature][badge-feature] Add the mean elements fitting and the epoch update to the
  two-body propagator: `fit_twobody_mean_elements`, `fit_twobody_mean_elements!`,
  `update_twobody_mean_elements_epoch`, `update_twobody_mean_elements_epoch!`, and the
  methods of `Propagators.fit_mean_elements` and `Propagators.fit_mean_elements!` for
  `Val(:TwoBody)`.
- ![Feature][badge-feature] The allocating functions that fit the mean elements and that
  update their epoch, `fit_*_mean_elements`, `update_*_mean_elements_epoch`, and
  `Propagators.fit_mean_elements`, accept the propagator constants with the same keyword of
  the initialization functions (`j2c`, `j4c`, and `m0`), whose number type selects the type
  of the fit, as `sgp4c` does in **SatelliteToolboxSgp4.jl**.
- ![Feature][badge-feature] The keyword `mean_elements_epoch` of every fitting function
  accepts a `DateTime` [UTC] besides a Julian Day, as `new_epoch` of the epoch update
  functions already did.
- ![Feature][badge-feature] Add the optional API functions `Propagators.propagator_data`,
  which returns the structure of the propagation theory wrapped by an `OrbitPropagator`,
  and `Propagators.is_initialized`, which tells whether a propagator created without
  initial elements has been initialized.
- ![Enhancement][badge-enhancement] The empty constructors of the propagator structures
  set the field `Δt` to `NaN`, which marks the structure as not initialized until an
  initialization function assigns its fields.
- ![Enhancement][badge-enhancement] Print the propagator structures `J2Propagator`,
  `J2OsculatingPropagator`, `J4Propagator`, `J4OsculatingPropagator`, and
  `TwoBodyPropagator` with the tree layout of **SatelliteToolboxBase.jl** v2.1, which
  follows the orbit data messages: the compact form shows the type and the epoch, and the
  rich form shows the sections with the initial mean elements and their epoch, the secular
  rates, the constants, and the last propagation instant. The `OrbitPropagator`
  wrappers print their type and name followed by the body of the wrapped structure, and an
  uninitialized propagator prints its status instead of undefined values.
- ![Enhancement][badge-enhancement] Replace **Crayons.jl** with **StyledStrings**, which
  only emits the terminal decorations when the output supports colors.
- ![Enhancement][badge-enhancement] Share the finite-difference and the ForwardDiff
  Jacobians, the epoch update, and the multi-threaded propagation of time vectors among the
  propagators, removing about 900 duplicated lines without changing the behavior or the
  allocation limits.
- ![Bugfix][badge-bugfix] Throw the documented `ArgumentError` when initializing a
  propagator with an invalid eccentricity. The conversion to the mean anomaly ran before the
  validation and raised a `DomainError` first.
- ![Bugfix][badge-bugfix] Fix the docstring examples of the fitting functions, which
  initialized the dummy propagator with an integer epoch and failed with an `InexactError`.
- ![Info][badge-info] Every function, including the private ones, now has a docstring, and
  the sources follow the coding style. The private Jacobian functions were merged into
  `_mean_elements_jacobian`.
- ![Info][badge-info] The propagator design description moved from `src/API.md` to the
  documentation, which now builds it as the page "API".
- ![Info][badge-info] Declare Aqua, JET, and AllocCheck as test dependencies instead of
  adding them to the test environment at run time. The quality checks moved to
  `test/quality.jl` and run on every stable Julia release.

Version 1.2.0
-------------

- ![Enhancement][badge-enhancement] The propagators now validate the eccentricity and the
  perigee radius during initialization, throwing an `ArgumentError` that points at the
  offending element instead of a `DomainError` raised by an internal square root.
- ![Enhancement][badge-enhancement] Improve the performance of the J2 and J4 osculating
  propagators by about 17%, and reduce the time to fit the mean elements with the
  finite-difference Jacobian by about half.
- ![Enhancement][badge-enhancement] Fitting the mean elements no longer allocates memory
  proportionally to the number of measurements when the position and velocity are passed as
  `Vector{Vector}`, which is the type used in all the documentation examples.
- ![Enhancement][badge-enhancement] The osculating Keplerian elements are now wrapped to
  [0, 2π), like the mean elements.
- ![Enhancement][badge-enhancement] The two-body propagator now wraps the mean anomaly to
  [0, 2π), improving the accuracy for large propagation times, especially when propagating
  in `Float32`.
- ![Enhancement][badge-enhancement] The J2 and J4 osculating propagators no longer compute
  a duplicated square root in the short-period corrections.
- ![Bugfix][badge-bugfix] Fix the short-period correction to the radial rate in the J2 and
  J4 osculating propagators, which used `(1 - e cos f)²` instead of `(1 + e cos f)²` in the
  term multiplying `sin(2u)`. The error vanishes for circular orbits and only affects the
  velocity.
- ![Bugfix][badge-bugfix] Fix `update_j2_mean_elements_epoch` and
  `update_j4_mean_elements_epoch`, which threw a `MethodError` when the elements were not
  `Float64`.
- ![Bugfix][badge-bugfix] Fix the `ntasks` keyword. It was rejected by the non-mutating
  vectorized functions, returned uninitialized elements when set to a value lower than one,
  and let surplus tasks write concurrently to the same output elements.
- ![Bugfix][badge-bugfix] Fix the epoch year encoding when obtaining the SGP4 mean elements,
  which returned an epoch one century away for the years before 1976.
- ![Bugfix][badge-bugfix] The functions that fit the mean elements now keep the epoch in
  `Tepoch` instead of converting it to the element type, and always return the documented
  `KeplerianElements{Tepoch, T}`.
- ![Bugfix][badge-bugfix] Fix the short-period corrections when the orbit crosses perigee,
  where a rounding difference could introduce a discontinuity, and avoid an overflow in the
  radial rate correction for orbits above roughly 50 900 km in `Float32`.
- ![Bugfix][badge-bugfix] Support arrays whose indices do not start at 1 when fitting the
  mean elements and when using the `OrbitStateVector` sinks.
- ![Bugfix][badge-bugfix] Fix the least-square fitting of the mean elements, which could
  never adjust a state component whose estimate was exactly zero, e.g. the z-axis
  components when the initial guess is an equatorial orbit. The correction limiter now
  bounds each state component against the norm of its position or velocity block instead
  of its own magnitude, which also removes a very slow convergence when a component is
  much smaller than the others, making the result robust across platforms.
- ![Info][badge-info] The short-period correction and the least-square algorithm that fits
  the mean elements are now implemented once and shared by the propagators, removing about
  1,100 duplicated lines.
- ![Info][badge-info] Fix many errors in the documentation, including wrong return types,
  wrong signatures, and the undocumented `jacobian_method` keyword.
- ![Info][badge-info] Fix the allocation tests of the Jacobian computed with
  `ForwardDiffJacobian`, which silently passed without testing anything because they used
  a stale keyword name.

Version 1.1.1
-------------

- ![Enhancement][badge-enhancement] Improve J2 and J4 propagator performance.
- ![Info][badge-info] Adopt BlueStyle and add repository instructions for coding agents.

Version 1.1.0
-------------

- ![Feature][badge-feature] The package now supports differentiability in all propagators.
  (PR [#5][gh-pr-5])

Version 1.0.0
-------------

- ![Info][badge-info] We dropped support for Julia 1.6. This version only supports the
  current Julia version and v1.10 (LTS).
- ![Info][badge-info] This version does not have breaking changes. We bump the version to
  1.0.0 because we now consider the API stable.

Version 0.3.3
-------------

- ![Feature][badge-feature] The functions `propagate`, `propagate!`, `propagate_to_epoch`,
  and `propagate_to_epoch!` of `Propagators` can now receive a sink option to change the
  type of the returned objects. We currently support `Tuple` (default) to keep the previous
  behavior or `OrbitStateVector` to return the results packed in an instance of
  `OrbitStateVector`.
- ![Feature][badge-feature] The functions to fit mean elements `fit_mean_elements` and
  `fit_mean_elements!` now supports `OrbitStateVector` as inputs.
- ![Deprecation][badge-deprecation] This package is no longer tested against Julia 1.6. The
  functions might work but there is no official support anymore. The next breaking release
  will remove the compatibility with Julia 1.6.

Version 0.3.2
-------------

- ![Feature][badge-feature] The functions `propagate`, `propagate!`, `propagate_to_epoch`,
  and `propagate_to_epoch!` of `Propagators` can now receive a vector of instants and the
  propagation will happen in multiple threads, if possible. The number of tasks can be
  configured using the keyword `ntasks`. If `ntasks = 1`, the algorithm falls back to the
  single thread version, leading to no overhead.
- ![Feature][badge-feature] We added support for the objects defined in `Dates` in the
  propagation functions. `propagate`, `propagate!`, and `step!` can now receive an object of
  type `Dates.Period` or `Dates.CompoundPeriod`. On the other hand, `propagate_to_epoch` and
  `propagate_to_epoch!` now supports an epoch specified using `DateTime`.
- ![Enhancement][badge-enhancement] We implemented dedicate `copy` to all propagators
  defined here, leading to a substantial gain compared to the previous version that relies
  on `deepcopy`.
- ![Enhancement][badge-enhancement] We increase the number of precompiled function
  signatures.

Version 0.3.1
-------------

- ![Enhancement][badge-enhancement] Minor source-code updates.
- ![Enhancement][badge-enhancement] We reduced the allocations in all functions that fit
  mean elements.

Version 0.3.0
-------------

- ![BREAKING][badge-breaking] We removed the possibility to add a mean motion
  time-derivative to J2 and J4 propagators. The theory we used to code those algorithms does
  not take into account such perturbations. Hence, the propagation accuracy would degrade
  very fast in those case with the mean motion perturbation.
- ![BREAKING][badge-breaking] The symbol to indicate a time-derivative in the structures of
  the propagators was changed from `δ` to `∂`.
- ![BREAKING][badge-breaking] ![Enhancement][badge-enhancement] We modified all the
  propagators to remove unnecessary variables in their structures after the redesign.
- ![Enhancement][badge-enhancement] We improved the J2 and J4 propagators given Kozai's
  theory.
- ![Enhancement][badge-enhancement] We updated the dependency compatibility bounds.

Version 0.2.1
-------------

- ![Bugfix][badge-bugfix] Fix the default constant name in the function `twobody`.

Version 0.2.0
-------------

- ![BREAKING][badge-breaking] ![Bugfix][badge-bugfix] The mean elements computed in the SGP4
  API `Propagators.mean_elements` was returning a set of osculating elements. This behavior
  is now fixed since the function now updates the epoch of the initial TLE to the last
  propagation instant.
- ![Feature][badge-feature] `Propagators` API now have functions to fit mean elements.
- ![Feature][badge-feature] We added functions to fit and update mean elements in all
  supported propagator but the two-body orbit propagator.
- ![Feature][badge-feature] We added the osculating version of the J4 propagator, called J4
  Osculation Orbit Propagator.
- ![Info][badge-info] The minimum version for **SatelliteToolboxSgp4.jl** was increased to
  2.1.

Version 0.1.0
-------------

- Initial version.
  - This version was based on the code in **SatelliteToolbox.jl**.

[badge-breaking]: https://img.shields.io/badge/Breaking-DC2626?style=flat-square
[badge-deprecation]: https://img.shields.io/badge/Deprecation-D97706?style=flat-square
[badge-feature]: https://img.shields.io/badge/Feature-16A34A?style=flat-square
[badge-enhancement]: https://img.shields.io/badge/Enhancement-0284C7?style=flat-square
[badge-bugfix]: https://img.shields.io/badge/Bugfix-DB2777?style=flat-square
[badge-info]: https://img.shields.io/badge/Info-475569?style=flat-square

[gh-pr-5]: https://github.com/JuliaSpace/SatelliteToolboxPropagators.jl/pull/5
