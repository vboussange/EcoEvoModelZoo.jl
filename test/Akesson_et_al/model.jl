using EcoEvoModelZoo
using ParametricModels
using Random; Random.seed!(1234)
using LinearAlgebra
using UnPack
using Test
using OrdinaryDiffEq
using Statistics
false ? (using PyPlot) : nothing

@testset "Akesson simple model" begin
    tstart = -4000 # starting time (relative to start of climate change at t = 0)
    tend = 2500 # time at which integration ends
    alg = Tsit5()
    abstol = 1e-6
    reltol = 1e-6
    # integrate prob with Julia
    tspan = (tstart, 0)
    saveat = tspan[1]:200:tspan[2]
    
    model = AkessonModel(ModelParams(;tspan, alg, reltol, abstol, saveat))

    sol = simulate(model)
    @test sol.retcode == :Success
end

@testset "Testing `eqs` integration against R script" begin
    R"source('/Users/victorboussange/ETHZ/projects/piecewise-inference/code/model/spatial_ecoevo-1.0.0/ecoevo_norun.R')"
    @rget SR SC S L rho kappa a eta eps W venv vmat s nmin aw bw Tmax Tmin Th arate Cmax Cmin tE d mig model ninit muinit
    ic = [ninit; muinit] # merge initial conditions into a vector
    S = S |> Int
    SR = SR |> Int
    SC = SC |> Int
    L = L |> Int
    x = collect(0:L-1) ./ (Float64(L) - 1.0)

    # coerce parameters into a dictionary
    pars = Dict{Symbol,Any}()
    @pack! pars = SR, SC, S, L, rho, kappa, a, eta, eps, W, venv, vmat, s, nmin, aw, bw, Tmax, Tmin, Th, arate, Cmax, Cmin, tE, d, mig, model, x
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

    uend = Array(before_cc)[:, :, end]

    # integrate prob with R
    println("Integrating with R")
    @time R"""
        before_cc <- ode(y=ic, times=seq(tstart, 0, by=200), func=eqs, parms=pars, method='bdf_d')
        uendR <- as.numeric(before_cc[nrow(before_cc),-1]) # final state -> new initial cond.
        """

    @rget uendR
    uend_rlike = [uend[1:SC+SR, :][:]; uend[SC+SR+1:end, :][:]]

    max_error_idx = argmax(abs.(uend_rlike .- uendR))

    @test abs((uend_rlike[max_error_idx] - uendR[max_error_idx]) / uendR[max_error_idx]) < 0.2

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