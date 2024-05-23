# Quick Start

```@meta
CurrentModule = SatelliteToolboxPropagators
```

```@setup quick_start
using SatelliteToolboxPropagators
```

Let's suppose we want to discover the ISS position on May 19, 2013, at 12:00:00. First, we
need to obtain a description of its orbit and use the corresponding propagator. Using the [Celestrak](https://celestrak.org) service, for example, we can obtain the ISS
[TLE](https://en.wikipedia.org/wiki/Two-line_element_set) on May 18, 2023:

```text
ISS (ZARYA)             
1 25544U 98067A   23138.86946505 -.00404279  00000+0 -74572-2 0  9994
2 25544  51.6431 113.8899 0006661 357.1286  88.0982 15.49924990397244
```

Since we have a TLE, we must use the SGP4/SDP4 propagator. Let's initialize the algorithm
first:

```@repl quick_start
iss_tle = tle"""
    ISS (ZARYA)
    1 25544U 98067A   23138.86946505 -.00404279  00000+0 -74572-2 0  9994
    2 25544  51.6431 113.8899 0006661 357.1286  88.0982 15.49924990397244"""

orbp = Propagators.init(Val(:SGP4), iss_tle)
```

Now we can propagate the ISS orbit to the desired epoch and obtain its state vector and mean
elements:

```@repl quick_start
Propagators.propagate_to_epoch!(orbp, date_to_jd(2023, 5, 19, 12, 0, 0))

Propagators.mean_elements(orbp)
```

!!! note

    Since we are using a TLE, the state vector is represented in the TEME reference frame.
