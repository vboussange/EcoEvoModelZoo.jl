#=
 This script defines ODE models using ParametricModels.jl
=#
using ParametricModels
using UnPack
using DocStringExtensions
using Statistics

"""
    $SIGNATURES

This model is inspired from [Huisman et al. 1999 Nature.](http://www.nature.com/articles/46540)
"""
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

This model is inspired from [Huisman et al. 1999 Nature.](http://www.nature.com/articles/46540), 
but where Leibig's law is replaced by 
imperfect substituable resources (smooth minimum). The smooth min function is parametrized by 
s, which is a trainable parameter.
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