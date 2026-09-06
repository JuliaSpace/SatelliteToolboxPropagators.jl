using Test

using Dates
using StaticArrays
using SatelliteToolboxPropagators

@testset "J2 Orbit Propagator" verbose = true begin
    include("./j2.jl")
end

@testset "J2 Osculating Orbit Propagator" verbose = true begin
    include("./j2osc.jl")
end

@testset "J4 Orbit Propagator" verbose = true begin
    include("./j4.jl")
end

@testset "J4 Osculating Orbit Propagator" verbose = true begin
    include("./j4osc.jl")
end

@testset "SGP4 Orbit Propagator" verbose = true begin
    include("./sgp4.jl")
end

@testset "Two-Body Orbit Propagator" verbose = true begin
    include("./twobody.jl")
end

@testset "Orbit Propagator API" verbose = true begin
    include("./api.jl")
end

# The quality and performance tools follow the Julia releases closely. Hence, we only run
# them on the stable releases.
if isempty(VERSION.prerelease)
    using Aqua
    using AllocCheck
    using ForwardDiff
    using JET

    @testset "Quality Tests" verbose = true begin
        include("./quality.jl")
    end

    @testset "Performance Tests" verbose = true begin
        include("./performance.jl")
    end
else
    @warn "The quality and performance tests are skipped on the Julia prereleases."
end
