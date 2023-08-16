SatelliteToolboxPropagator.jl Changelog
=======================================

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
