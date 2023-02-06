using DocStringExtensions
using Distributions
using LinearAlgebra
"""
    $SIGNATURES

Returns a Dictionary with all params for Akesson model
"""
function init_params(;
    S=50,
    L=50, # number of patches
    vbar=0.1,
    dbar=1e-5,
    model_type="Tdep_trophic",
    kappa=0.1, # intrinsic mortality parameter
    eta=1, # competition width (centigrade; only for Tdep and Tdep_trophic)
    nmin=1e-5, # below this threshold density, genetic variances are reduced
    aw=0.1, # (negative) slope of trait-dependence of tolerance width
    bw=4, # intercept of trait-dependence of tolerance width
    Tmax=25.0, # initial mean temperature at equator
    Tmin=-10.0,# initial mean temperature at poles
    Cmax=9.66, # projected temperature increase at poles
    Cmin=1.26 # projected temperature increase at equator
)
    SR = S # number of resource species

    if (model_type in ["trophic", "Tdep_trophic"])
        SC = S # ...consumer species
    else
        SC = 0 # number of consumer species: 0, unless we have...
    end
    S = SR + SC # set S to be the total number of species

    # random- and trophic-dependent quantities
    v = rand(Uniform(0.5 * vbar, 1.5 * vbar), SR) # resource genetic variances

    d = rand(Uniform(0.1 * dbar, 10.0 * dbar), SR) # resource dispersal rates
    rho = rand(Uniform(0.9, 1.1), SR)  # resource growth-tolerance tradeoff parameter
    a = zeros(S, S) # initialize full competition matrix (resources+consumers)
    aP = rand(Uniform(0.15 * 0.5, 0.15 * 1.5), SR, SR) # resource comp coeffs
    aP[diagind(aP)] .= rand(Uniform(0.2 * 0.5, 0.2 * 1.5), SR)  # resource intraspecific comp coeffs
    a[1:SR, 1:SR] = aP # top left block: resources
    W = zeros(S, S) # create feeding network: nothing if no consumers
    Th = ones(S) # handling times in type II f.r. (dummy value if no consumers)
    arate = ones(S) # attack rates in type II f.r. (dummy value if no

    if model_type in ["trophic", "Tdep_trophic"]
        v = [v; rand(Uniform(0.5 * vbar, 1.5 * vbar), SC)] # add consumer genetic variances
        d = [d; rand(Uniform(0.1 * dbar, 10.0 * dbar), SC)] # add consumer dispersal rates
        rho = [rho; rand(Uniform(0.9*0.1, 1.1*0.1), SC)] # add consumer tradeoff parameters
        aH = zeros(SC, SC) # initialize competition matrix (consumers)
        a[SR+1:S, SR+1:S] .= aH # bottom right: consumers
        W = generate_network(SR, SC) # trophic feeding network
        Th[SR+1:S] .= rand(SC) .* (1 .- 0.5) .+ 0.5 # handling times in type II f.r.
        arate[SR+1:S] = rand(SC) .* (10 .- 1) .+ 1 # attack rates in type II f.r.
    end

    # all other parameters
    venv = vbar # environmental variance
    vmat = repeat(v, outer=(1, L)) # genetic variances at each patch
    s = v .+ venv # species' total phenotypic variances
    eps = vcat(zeros(SR), 0.3 * ones(SC)) # feeding efficiency of consumers
    tE = 300.0 # time at which climate change stops (assuming it starts at t = 0)

    # dispersal matrix
    mig = zeros(L, L) # initialize dispersal matrix
    for k in 2:L
        mig[k-1, k] = 1 # each species can only migrate to the two
    end
    mig = mig + mig' # nearest-neighbor patches

    # initial temperatures
    if L > 1
        Tempinit = Temp.(range(0, 1, L), 0, tE, Cmax, Cmin, Tmax, Tmin)
        x = collect(0:L-1) ./ (Float64(L) - 1.0)
    else
        Tempinit = [0.5]
        x = [0.5]
    end


    # initial conditions
    muinit = [Tmin + (Tmax - Tmin) * i / SR for i in 1:SR, j in 1:L] # initial trait means
    ninit = zeros(S, L) # reserve memory for initial densities
    for i in 1:SR
        ninit[i, :] = exp.(-(muinit[i, 1] .- Tempinit) .^ 2 / (2 * 2^2))
    end

    # initial traits and densities for consumers
    if model_type in ["trophic", "Tdep_trophic"]
        muinit = vcat(muinit, [Tmin + (Tmax - Tmin) * i / SC for i in 1:SC, j in 1:L])
        for i in (SR+1):S
            ninit[i, :] = exp.(-(muinit[i, 1] .- Tempinit) .^ 2 / (2 * 2^2))
        end
    end

    ic = [ninit; muinit] # merge initial conditions into a vector


    # coerce parameters into a dictionary
    pars = Dict{Symbol,Any}()
    @pack! pars = SR, SC, S, L, rho, kappa, a, eta, eps, W, venv, vmat, s, nmin, aw, bw, Tmax, Tmin, Th, arate, Cmax, Cmin, tE, d, mig, model_type, x

    return pars, ic
end