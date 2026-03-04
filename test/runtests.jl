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

if isempty(VERSION.prerelease)
    using Pkg

    Pkg.add("JET")
    Pkg.add("AllocCheck")
    Pkg.add("Aqua")
    Pkg.add("Lux")
    Pkg.add("Optimisers")
    Pkg.add("Zygote")

    using JET
    using AllocCheck
    using Aqua

    @testset "Performance Tests" verbose = true begin
        include("./performance.jl")
    end

    @testset "Lux Extension (ML-∂SGP4)" verbose = true begin
        include("./extension.jl")
    end
else
    @warn "Performance checks and extensions not guaranteed to work on julia-nightly, skipping tests"
end
