SatelliteToolboxPropagator.jl Changelog
=======================================

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

[badge-breaking]: https://img.shields.io/badge/BREAKING-red.svg
[badge-deprecation]: https://img.shields.io/badge/Deprecation-orange.svg
[badge-feature]: https://img.shields.io/badge/Feature-green.svg
[badge-enhancement]: https://img.shields.io/badge/Enhancement-blue.svg
[badge-bugfix]: https://img.shields.io/badge/Bugfix-purple.svg
[badge-info]: https://img.shields.io/badge/Info-gray.svg
