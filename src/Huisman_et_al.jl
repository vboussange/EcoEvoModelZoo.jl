#=
 This script defines ODE models using ParametricModels.jl
=#
using ParametricModels
using UnPack
using DocStringExtensions
using Statistics

"""
    $SIGNATURES

The `ResourceCompetition` model is a type of ecological model that simulates
competition for resources between different species of plankton. This model is
inspired by Huisman et al. 1999 Nature.

# Arguments
- `mp` is a `ModelParams` object that contains information about the time span, the algorithm to use for numerical integration, and the tolerances to use.
- `nN` is the number of plankton species.
- `nR` is the number of resources.
- `S` is a vector of length nR containing the supply rates of each resource.
- `mu`: function for the ResourceCompetition model is set to the default function:

```julia
μ(R, r, K) = r .* R ./ (K .+ R)
```
where:

- `R` is a vector of length nR containing the current concentration of each resource.
- `r` is a vector of length nN containing the maximum uptake rate of each plankton species.
- `K` is a vector of length nN containing the half-saturation constant of each plankton species."""
Base.@kwdef struct ResourceCompetition{MP,NN,NR,Mu,SS} <: AbstractModel
    mp::MP
    nN::NN # number of planktons
    nR::NR # number of resources
    S::SS # supply rates
    mu::Mu # mu function
end

ResourceCompetition(mp, nN, nR, S) = ResourceCompetition(mp, nN, nR, S, μ)

function (model::ResourceCompetition)(du, u, p, t)
    @unpack r, K, C, m, D = p
    @unpack nN, nR, mu, S = model
    ũ = max.(u, 0.)
    dN = @view du[1:nN]; dR = @view du[nN + 1:end]
    N = @view ũ[1:nN]; R = @view ũ[nN + 1:end]

    _mu = mu(R, r, K)
    # @assert !any(isinf.(_mu) .|| isnan.(_mu)) string("problem with mu", @show ForwardDiff.value(mu(R, r, K)), @show ForwardDiff.value(R))

    dN .= N .* (_mu .- m[])
    dR .= D[] .* (S .- R) .- sum(transpose(C .* _mu .* N), dims=2)

    # @assert !any(isinf.(dN) .|| isnan.(dN)) string("problem with dN")
    # @assert !any(isinf.(dR) .|| isnan.(dR)) string("problem with dR")

end

"""
    $SIGNATURES

This model is inspired from [Huisman et al. 1999
Nature.](http://www.nature.com/articles/46540), but where Leibig's law is
replaced by imperfect substituable resources ([smooth
minimum](https://en.wikipedia.org/wiki/Smooth_maximum)). The smooth min function
is parametrized by s, which is a trainable parameter.

# Arguments
- `mp::MP`: ModelParams struct containing parameter values for the model.
- `nN::NN`: Number of plankton species in the model.
- `nR::NR`: Number of resource species in the model.
- `S::SS`: Supply rates.

# Example
```
alg = BS3()
abstol = 1e-6
reltol = 1e-6
tspan = (0.,1000.)
step = 2.
nN = 5
nR = 5
r = ones(nN)
m = D = 0.25
S = [6., 10., 14., 4, 9]
K = [0.39 0.34 0.30 0.24 0.23;
    0.22 0.39 0.34 0.30 0.27;
    0.27 0.22 0.39 0.34 0.30;
    0.30 0.24 0.22 0.39 0.34;
    0.34 0.30 0.22 0.20 0.39]'

C = [0.04 0.04 0.07 0.04 0.04;
    0.08 0.08 0.08 0.10 0.08;
    0.10 0.10 0.10 0.10 0.14;
    0.05 0.03 0.03 0.03 0.03;
    0.07 0.09 0.07 0.07 0.07]'
s = -7.

p_true = (r = r, m = [m], D = [D], K = K, C = C, s = Float64[s])

N0 = [0.1 + i / 100 for i in 1:5]
R0 = S
u0 = [N0;R0]
tsteps = tspan[1]:step:tspan[2]

model = ResourceCompetitionSmoothMin(ModelParams(;p = p_true,
                                        tspan,
                                        u0,
                                        alg,
                                        reltol,
                                        abstol,
                                        saveat = tsteps
                                        ), 
                                        nN, nR, S)

# Species 1:5
sol = simulate(model, u0 = u0)
```
"""
Base.@kwdef struct ResourceCompetitionSmoothMin{MP,NN,NR,SS} <: AbstractModel
    mp::MP
    nN::NN # number of planktons
    nR::NR # number of resources
    S::SS # supply rates
end

function (model::ResourceCompetitionSmoothMin)(du, u, p, t)
    @unpack r, K, C, m, D, s = p
    @unpack nN, nR, S = model
    ũ = max.(u, 0.)
    dN = @view du[1:nN]; dR = @view du[nN + 1:end]
    N = @view ũ[1:nN]; R = @view ũ[nN + 1:end]

    _mu = μsmooth(R, r, K, s)
    # @assert !any(isinf.(_mu) .|| isnan.(_mu)) string("problem with mu", @show ForwardDiff.value(mu(R, r, K)), @show ForwardDiff.value(R))

    dN .= N .* (_mu .- m[])
    dR .= D[] .* (S .- R) .- sum(transpose(C .* _mu .* N), dims=2)

    # @assert !any(isinf.(dN) .|| isnan.(dN)) string("problem with dN")
    # @assert !any(isinf.(dR) .|| isnan.(dR)) string("problem with dR")

end

μ(R, r, K) = minimum(r .* R' ./ (K .+ R'), dims = 2)
smoothmin(s, X) = sum(X .* exp.(s .* X), dims=2) ./ sum(exp.(s .* X), dims=2)
μsmooth(R, r, K, s) = smoothmin(s, r .* R' ./ (K .+ R'))