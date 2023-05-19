Quick Start
===========

```@meta
CurrentModule = SatelliteToolboxPropagators
DocTestSetup = quote
    using SatelliteToolboxPropagators
end
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

```jldoctest SGP4
julia> iss_tle = tle"""
       ISS (ZARYA)             
       1 25544U 98067A   23138.86946505 -.00404279  00000+0 -74572-2 0  9994
       2 25544  51.6431 113.8899 0006661 357.1286  88.0982 15.49924990397244"""
TLE:
                     Name : ISS (ZARYA)
         Satellite number : 25544
 International designator : 98067A
       Epoch (Year / Day) : 23 / 138.86946505 (2023-05-18T20:52:01.780)
       Element set number : 999
             Eccentricity :   0.00066610 deg
              Inclination :  51.64310000 deg
                     RAAN : 113.88990000 deg
      Argument of perigee : 357.12860000 deg
             Mean anomaly :  88.09820000 deg
          Mean motion (n) :  15.49924990 revs/day
        Revolution number : 39724
                       B* : -0.007457 1/[er]
                    ṅ / 2 : -0.004043 rev/day²
                    n̈ / 6 : 0.000000 rev/day³

julia> orbp = Propagators.init(Val(:SGP4), iss_tle)
OrbitPropagatorSgp4{Float64, Float64}:
   Propagator name : SGP4 Orbit Propagator
  Propagator epoch : 2023-05-18T20:52:01.780
  Last propagation : 2023-05-18T20:52:01.780
```

Now we can propagate the ISS orbit to the desired epoch and obtain its state vector and mean
elements:

```jldoctest SGP4
julia> Propagators.propagate_to_epoch!(orbp, date_to_jd(2023, 5, 19, 12, 0, 0))
([-2.7473591812349055e6, 6.194023092402147e6, 470636.34213714604], [-4186.242303388817, -2319.0180933770293, 5989.510493513584])

julia> Propagators.mean_elements(orbp)
KeplerianElements{Float64, Float64}:
           Epoch :   2.46008e6 (2023-05-18T21:07:09.751)
 Semi-major axis :   6.79714     km
    Eccentricity :   0.000670975
     Inclination :  51.6627      °
            RAAN : 110.771       °
 Arg. of Perigee : 359.796       °
    True Anomaly :   5.97883     °
```

!!! note
    Since we are using a TLE, the state vector is represented in the TEME reference frame.
