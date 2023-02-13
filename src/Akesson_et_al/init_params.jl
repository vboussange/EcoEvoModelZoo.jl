using DocStringExtensions
using Distributions
using LinearAlgebra
"""
    $SIGNATURES

Returns a Dictionary with all *biological* parameters for Akesson model.

# Arguments
- `model_type`: can take values "normal", "Tdep", "Tdep_trophic" and "trophic".
"""
function init_params_akesson_model(land::Landscape, 
                            temp::Temperature,
                            SR, SC;
                            width_growth::WG = WidthGrowth{:TraitDep}(), 
                            competition::CP = Competition{:TraitDep}(), 
                            trophic::Tr = Trophic{false}(),
                            vbar=0.1, # mean genetic variance, used to generate genetic variance vector
                            dbar=1e-5, # mean dispersal rate, used to generate dispersal vector
                            kappa=[0.1], # intrinsic mortality parameter
                            eta=[1.], # competition width (centigrade; only for Tdep and Tdep_trophic)
                            nmin=1e-5, # below this threshold density, genetic variances are reduced
                            aw=[0.1], # (negative) slope of trait-dependence of tolerance width
                            bw=[4.], # intercept of trait-dependence of tolerance width
                            w = [0.2], # tolerance width for when width_growth = WidthGrowth{:Standard}()
                        ) where {WG, CP, Tr}

    @unpack L, x = land
    # coerce parameters into a dictionary
    pars = Dict{Symbol,Any}()
    S = SR + SC # set S to be the total number of species

    # random- and trophic-dependent quantities
    v = rand(Uniform(0.5 * vbar, 1.5 * vbar), SR) # resource genetic variances
    d = rand(Uniform(0.1 * dbar, 10.0 * dbar), SR) # resource dispersal rates
    rho = rand(Uniform(0.9, 1.1), SR)  # resource growth-tolerance tradeoff parameter

    if (CP <: Competition{:Standard})
        a = zeros(S, S) # initialize full competition matrix (resources+consumers)
        aP = rand(Uniform(0.15 * 0.5, 0.15 * 1.5), SR, SR) # resource comp coeffs
        aP[diagind(aP)] .= rand(Uniform(0.2 * 0.5, 0.2 * 1.5), SR)  # resource intraspecific comp coeffs
        a[1:SR, 1:SR] = aP # top left block: resources
        if Tr <: Trophic{true}
            aH = zeros(SC, SC) # initialize competition matrix (consumers)
            a[SR+1:S, SR+1:S] .= aH # bottom right: consumers
        end
        @pack! pars = a
    elseif (CP <: Competition{:TraitDep})
        @pack! pars = eta
    end


    if Tr <: Trophic{true}
        v = [v; rand(Uniform(0.5 * vbar, 1.5 * vbar), SC)] # add consumer genetic variances
        d = [d; rand(Uniform(0.1 * dbar, 10.0 * dbar), SC)] # add consumer dispersal rates
        rho = [rho; rand(Uniform(0.9*0.1, 1.1*0.1), SC)] # add consumer tradeoff parameters
        W = generate_network(SR, SC) # trophic feeding network
        Th = ones(S) # handling times in type II f.r. (dummy value if no consumers)
        arate = ones(S) # attack rates in type II f.r. (dummy value if no
        Th[SR+1:S] .= rand(SC) .* (1 .- 0.5) .+ 0.5 # handling times in type II f.r.
        arate[SR+1:S] = rand(SC) .* (10 .- 1) .+ 1 # attack rates in type II f.r.
        eps = vcat(zeros(SR), 0.3 * ones(SC)) # feeding efficiency of consumers
        @pack! pars = W, arate, eps
    end

    if WG <: WidthGrowth{:TraitDep}
        @pack! pars = aw, bw
    elseif WG <: WidthGrowth{:Standard}
        @pack! pars = w
    end
    # all other parameters
    venv = [vbar] # environmental variance

    # initial conditions
    # initial temperatures
    if L > 1 
        Tempinit = temp.(x, 0)
    else
        Tempinit = [0.5]
    end

    @unpack Tmax, Tmin = temp
    muinit = [Tmin + (Tmax - Tmin) * i / SR for i in 1:SR, j in 1:L] # initial trait means
    ninit = zeros(S, L) # reserve memory for initial densities
    for i in 1:SR
        ninit[i, :] = exp.(-(muinit[i, 1] .- Tempinit) .^ 2 / (2 * 2^2))
    end

    # initial traits and densities for consumers
    if Tr <: Trophic{true}
        muinit = vcat(muinit, [Tmin + (Tmax - Tmin) * i / SC for i in 1:SC, j in 1:L])
        for i in (SR+1):S
            ninit[i, :] = exp.(-(muinit[i, 1] .- Tempinit) .^ 2 / (2 * 2^2))
        end
    end

    ic = [ninit; muinit] # merge initial conditions into a vector
    
    @pack! pars = rho, kappa, venv, v, nmin, d

    return pars, ic, Tempinit, x
end