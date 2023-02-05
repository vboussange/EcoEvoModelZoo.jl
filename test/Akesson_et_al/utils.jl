cd(@__DIR__)
using Test
using UnPack
using RCall
using Random
using Distributions
using LinearAlgebra
using ForwardDiff
using OrdinaryDiffEq
using EcoEvoModelZoo


@testset "generate_network" begin
    SR = 10
    SC = 10
    @test size(generate_network(SR, SC)) == (SR + SC, SR + SC)
end

@testset "organize_data" begin
    println("no test implemented")
end

@testset "smoothstep" begin
    @test smoothstep(-1.0) == 0
    @test smoothstep(2.0) == 1
    @test 0 < smoothstep(0.5) < 1
end

@testset "Temp" begin
    Tmax = 25.0 # initial mean temperature at equator
    Tmin = -10.0 # initial mean temperature at poles
    Cmax = 9.66 # projected temperature increase at poles
    Cmin = 1.26 # projected temperature increase at equator
    x = 0.5
    t = 50
    tE = 300.0 # time at which climate change stops (assuming it starts at t = 0)
    @test Temp(x, t, tE, Cmax, Cmin, Tmax, Tmin) isa Number
    x = 1:10
    @test Temp.(x, t, tE, Cmax, Cmin, Tmax, Tmin) isa Vector

end

@testset "funcresp" begin
    SR = 10
    SC = 10
    n = rand(SC + SR) # Vector of population densities of all species in a given patch
    W = generate_network(SR, SC)
    Th = ones(SC + SR) # Vector of handling times (with dummy values for resource species)
    arate = ones(SC + SR)
    F = zeros(SR + SC, SR + SC)
    funcresp!(F, n, Th, arate, W)
    @test F isa Matrix
    @test !any(isnan.(F))
    @test !any(isinf.(F))
end

@testset "init_params" begin

    # default parameters
    pars, ic = init_params()

    # default R parameters
    R"source('/Users/victorboussange/ETHZ/projects/piecewise-inference/code/model/spatial_ecoevo-1.0.0/ecoevo_norun.R')"
    @rget SR SC S L rho kappa a eta eps W venv vmat s nmin aw bw Tmax Tmin Th arate Cmax Cmin tE d mig model ninit muinit
    ic = [ninit; muinit] # merge initial conditions into a vector
    S = S |> Int
    SR = SR |> Int
    SC = SC |> Int
    L = L |> Int
    x = collect(0:L-1) ./ (Float64(L) - 1.0)

    # coerce parameters into a dictionary
    parsR = Dict{String,Any}()
    @pack! parsR = SR, SC, S, L, rho, kappa, a, eta, eps, W, venv, vmat, s, nmin, aw, bw, Tmax, Tmin, Th, arate, Cmax, Cmin, tE, d, mig, model, x

    @test pars[:S] == parsR["S"]
    @test pars[:SR] == parsR["SR"]
    @test pars[:SC] == parsR["SC"]
    @test pars[:L] == parsR["L"]
    @test size(pars[:rho]) == size(parsR["rho"])
    # uncomplete
end

@testset "Basic tests of `eqs`" begin

    tstart = -4000 # starting time (relative to start of climate change at t = 0)
    tE = 300.0 # time at which climate change stops (assuming it starts at t = 0)
    tend = 2500 # time at which integration ends
    pars, ic = init_params()

    pars = NamedTuple([pair for pair in pars])

    dudt = similar(ic)

    eqs!(dudt, ic, pars, 0.0)

    @test !any(isnan.(dudt))
    @test !any(isinf.(dudt))
end


@testset "Compare `eqs` evaluation against R script" begin
    R"source('/Users/victorboussange/ETHZ/projects/piecewise-inference/code/model/spatial_ecoevo-1.0.0/ecoevo_norun.R')"
    @rget SR SC S L rho kappa a eta eps W venv vmat s nmin aw bw Tmax Tmin Th arate Cmax Cmin tE d mig model ninit muinit
    ic = [ninit; muinit] # merge initial conditions into a vector
    S = S |> Int
    SR = SR |> Int
    SC = SC |> Int
    L = L |> Int
    x = collect(0:L-1) ./ (Float64(L) - 1.0)

    # coerce parameters into a dictionary
    pars = Dict{String,Any}()
    @pack! pars = SR, SC, S, L, rho, kappa, a, eta, eps, W, venv, vmat, s, nmin, aw, bw, Tmax, Tmin, Th, arate, Cmax, Cmin, tE, d, mig, model, x

    dudt = similar(ic)

    eqs!(dudt, ic, pars, 10.0)

    # for comparing, one needs to reorder
    dudt_rlike = [dudt[1:SC+SR, :][:]; dudt[SC+SR+1:end, :][:]]

    # test with R
    dudtR = R"eqs(1., ic, pars)"[1]

    @test all(isapprox.(dudt_rlike, dudtR, atol=1e-3))
end

@testset "eqs with L = 1" begin

    tstart = -4000 # starting time (relative to start of climate change at t = 0)
    tend = 2500 # time at which integration ends
    pars, ic = init_params()
    pars = NamedTuple([pair for pair in pars])

    dudt = similar(ic)

    eqs!(dudt, ic, pars, 10.0)

    @test !any(isnan.(dudt))
    @test !any(isinf.(dudt))

end

@testset "Testing with L = 1 (no spatial structure)" begin

    tstart = -4000 # starting time (relative to start of climate change at t = 0)
    tend = 2500 # time at which integration ends
    pars, ic = init_params(L=1, model="normal")
    pars = NamedTuple([pair for pair in pars])


    # integrate prob with Julia
    tspan = (tstart, 0)
    saveat = tspan[1]:200:tspan[2]
    prob = ODEProblem(eqs!, ic, tspan, pars)

    println("Integrating with Julia")
    @time before_cc = solve(prob,
        alg=Tsit5(),
        saveat=saveat)

    @test before_cc.retcode == :Success

end

@testset "Testing differentiation, in the setting with no spatial structure" begin
    tstart = -4000 # starting time (relative to start of climate change at t = 0)
    tend = 2500 # time at which integration ends
    pars, ic = init_params(L=1, S = 3, model="normal")

    @unpack SR, SC, S, L, rho, kappa, a, eta, eps, W, venv, vmat, s, nmin, aw, bw, Tmax, Tmin, Th, arate, Cmax, Cmin, tE, d, mig, model, x = pars

    function dummy_loss(a)
        # coerce parameters into a dictionary
        pars = Dict{Symbol,Any}()
        @pack! pars = SR, SC, S, L, rho, kappa, a, eta, eps, W, venv, vmat, s, nmin, aw, bw, Tmax, Tmin, Th, arate, Cmax, Cmin, tE, d, mig, model, x

        pars = NamedTuple([pair for pair in pars])
        prob = ODEProblem(eqs!, ic, tspan, pars)

        before_cc = solve(prob,
                        alg=Tsit5(),
                        saveat=saveat)
       sum(before_cc.^2)
    end

    # Checking differentiation
    J = ForwardDiff.gradient(dummy_loss, a)

    @test F isa Matrix
    @test !any(isnan.(F))
    @test !any(isinf.(F))
end