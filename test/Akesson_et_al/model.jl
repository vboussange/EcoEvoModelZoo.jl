using EcoEvoModelZoo
using ParametricModels
using Random; Random.seed!(1234)
using LinearAlgebra
using UnPack
using Test
using OrdinaryDiffEq
using Statistics
using ForwardDiff

test_rcall = true
test_plot = false

test_rcall ? (using RCall) : nothing
test_plot ? (using PyPlot) : nothing

@testset "Init model" begin
    model = init_akesson_model()
end

@testset "Testing Akesson default model" begin
    tstart = -4000 # starting time (relative to start of climate change at t = 0)
    tend = 2500 # time at which integration ends
    alg = Tsit5()
    abstol = 1e-6
    reltol = 1e-6
    # integrate prob with Julia
    tspan = (tstart, 0)
    saveat = tspan[1]:200:tspan[2]
    
    model = init_akesson_model(mp = ModelParams(;tspan, alg, reltol, abstol, saveat))

    sol = simulate(model)
    @test SciMLBase.successful_retcode(sol.retcode)
end

if test_rcall
    @testset "Testing `eqs` integration against R script - Tdep_trophic" begin
        R"source('/Users/victorboussange/ETHZ/projects/piecewise-inference/code/model/spatial_ecoevo-1.0.0/ecoevo_norun.R')"
        @rget SR SC S L rho kappa a eta eps W venv v s nmin aw bw Tmax Tmin Th arate Cmax Cmin tE d mig model ninit muinit
        ic = [ninit; muinit] # merge initial conditions into a vector
        S = S |> Int
        SR = SR |> Int
        SC = SC |> Int
        L = L |> Int
        x = collect(0:L-1) ./ (Float64(L) - 1.0)

        # integrate prob with Julia
        alg = Tsit5()
        abstol = 1e-6
        reltol = 1e-6
        tstart = -4000 # starting time (relative to start of climate change at t = 0)
        tspan = (tstart, 0)
        saveat = tspan[1]:200:tspan[2]
        
         # coerce parameters into a dictionary
        pars = Dict{Symbol,Any}()
        @pack! pars = rho, kappa, a, eta, eps, W, venv, v, nmin, aw, bw, Th, arate, d
        pars = NamedTuple([pair for pair in pars])
        model = init_akesson_model(;SR,
                                    SC,                         
                                    mp = ModelParams(;p = pars, u0 = ic, tspan, alg, reltol, abstol, saveat),
                                    land = Landscape(L, x, mig),
                                    temp = Temperature(;Tmax, Tmin, Cmax, Cmin, tE,),
                                    width_growth = WidthGrowth{:TraitDep}(), 
                                    competition = Competition{:TraitDep}(), 
                                    trophic= Trophic{true}(), 
                                    )

        println("Integrating with Julia")
        @time before_cc = simulate(model)
        @test SciMLBase.successful_retcode(before_cc.retcode)

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
    @testset "Testing `eqs` integration against R script - Tdep" begin
        R"source('/Users/victorboussange/ETHZ/projects/EcoEvoModelZoo.jl/test/Akesson_et_al/ecoevo_norun_Tdep.R')"
        @rget SR SC S L rho kappa a eta eps W venv v nmin aw bw Tmax Tmin Th arate Cmax Cmin tE d mig model ninit muinit
        ic = [ninit; muinit] # merge initial conditions into a vector
        S = S |> Int
        SR = SR |> Int
        SC = SC |> Int
        L = L |> Int
        x = collect(0:L-1) ./ (Float64(L) - 1.0)

        # integrate prob with Julia
        alg = Tsit5()
        abstol = 1e-6
        reltol = 1e-6
        tstart = -4000 # starting time (relative to start of climate change at t = 0)
        tspan = (tstart, 0)
        saveat = tspan[1]:200:tspan[2]
        
         # coerce parameters into a dictionary
        pars = Dict{Symbol,Any}()
        @pack! pars = rho, kappa, a, eta, eps, W, venv, v, nmin, aw, bw, Th, arate, d
        pars = NamedTuple([pair for pair in pars])
        model = init_akesson_model(;SR,
                                    SC,                         
                                    mp = ModelParams(;p = pars, u0 = ic, tspan, alg, reltol, abstol, saveat),
                                    land = Landscape(L, x, mig),
                                    temp = Temperature(;Tmax, Tmin, Cmax, Cmin, tE,),
                                    width_growth = WidthGrowth{:TraitDep}(), 
                                    competition = Competition{:TraitDep}(), 
                                    trophic= Trophic{false}(), 
                                    )

        println("Integrating with Julia")
        @time before_cc = simulate(model)
        @test SciMLBase.successful_retcode(before_cc.retcode)

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
    @testset "Testing `eqs` integration against R script - standard model" begin
        R"source('/Users/victorboussange/ETHZ/projects/EcoEvoModelZoo.jl/test/Akesson_et_al/ecoevo_norun_std.R')"
        @rget SR SC S L rho kappa a eta eps W venv v nmin aw bw Tmax Tmin Th arate Cmax Cmin tE d mig model ninit muinit
        ic = [ninit; muinit] # merge initial conditions into a vector
        S = S |> Int
        SR = SR |> Int
        SC = SC |> Int
        L = L |> Int
        x = collect(0:L-1) ./ (Float64(L) - 1.0)

        # integrate prob with Julia
        alg = Tsit5()
        abstol = 1e-6
        reltol = 1e-6
        tstart = -4000 # starting time (relative to start of climate change at t = 0)
        tspan = (tstart, 0)
        saveat = tspan[1]:200:tspan[2]
        
         # coerce parameters into a dictionary
        pars = Dict{Symbol,Any}()
        V = s
        @pack! pars = rho, kappa, a, eta, eps, W, venv, v, nmin, aw, bw, Th, arate, d
        pars = NamedTuple([pair for pair in pars])
        model = init_akesson_model(;SR,
                                    SC,                         
                                    mp = ModelParams(;p = pars, u0 = ic, tspan, alg, reltol, abstol, saveat),
                                    land = Landscape(L, x, mig),
                                    temp = Temperature(;Tmax, Tmin, Cmax, Cmin, tE,),
                                    width_growth = WidthGrowth{:TraitDep}(), 
                                    competition = Competition{:Standard}(), 
                                    trophic= Trophic{false}(), 
                                    )

        println("Integrating with Julia")
        @time before_cc = simulate(model)
        @test SciMLBase.successful_retcode(before_cc.retcode)

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
end

@testset "Testing with L = 1 (no spatial structure)" begin
    tstart = -4000 # starting time (relative to start of climate change at t = 0)
    tend = 2500 # time at which integration ends
    alg = Tsit5()
    abstol = 1e-6
    reltol = 1e-6
    # integrate prob with Julia
    tspan = (tstart, 0)
    saveat = tspan[1]:200:tspan[2]
    
    model = init_akesson_model(mp = ModelParams(;tspan, alg, reltol, abstol, saveat);)

    sol = simulate(model)
    @test SciMLBase.successful_retcode(sol.retcode)
end

@testset "Testing differentiation, in the setting with no spatial structure" begin

    tstart = -4000 # starting time (relative to start of climate change at t = 0)
    tend = 2500 # time at which integration ends
    alg = Tsit5()
    abstol = 1e-6
    reltol = 1e-6
    # integrate prob with Julia
    tspan = (tstart, 0)
    saveat = tspan[1]:200:tspan[2]
    
    model = init_akesson_model(mp = ModelParams(;tspan, alg, reltol, abstol, saveat);
                                land = Landscape(1),
                                temp = Temperature(),
                                width_growth = WidthGrowth{:Standard}(), 
                                competition = Competition{:Standard}(), 
                                trophic= Trophic{false}(), )

    function dummy_loss(rho)
        before_cc = simulate(model, p = (rho = rho,))
        sum(before_cc.^2)
    end

    # Checking differentiation
    rho = get_p(model).rho
    J = ForwardDiff.gradient(dummy_loss, rho)

    @test J isa Vector
    @test !any(isnan.(J))
    @test !any(isinf.(J))
end

@testset "Testing differentiation of eta" begin

    tstart = -4000 # starting time (relative to start of climate change at t = 0)
    tend = 2500 # time at which integration ends
    alg = Tsit5()
    abstol = 1e-6
    reltol = 1e-6
    # integrate prob with Julia
    tspan = (tstart, 0)
    saveat = tspan[1]:200:tspan[2]
    
    model = init_akesson_model(mp = ModelParams(;tspan, alg, reltol, abstol, saveat);
                                land = Landscape(1),
                                temp = Temperature(),
                                width_growth = WidthGrowth{:TraitDep}(), 
                                competition = Competition{:TraitDep}(), 
                                trophic= Trophic{false}(), )

    function dummy_loss(eta)
        before_cc = simulate(model, p = (eta = eta,))
        sum(before_cc.^2)
    end

    # Checking differentiation
    rho = get_p(model).eta
    J = ForwardDiff.gradient(dummy_loss, eta)

    @test J isa Vector
    @test !any(isnan.(J))
    @test !any(isinf.(J))
end