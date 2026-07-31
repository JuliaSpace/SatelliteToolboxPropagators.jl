SatelliteToolboxPropagator.jl Changelog
=======================================

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
  components when the initial guess is an equatorial orbit.
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
