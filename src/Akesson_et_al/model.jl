using ParametricModels
using UnPack
using DocStringExtensions

"""
    $SIGNATURES

This model is inspired from [Akesson et al. 2021](https://www.nature.com/articles/s41467-021-24977-x).
General version of the model, corresponding to original R script.

Specialized variants are provided in other `AbstractModel`s.
"""
Base.@kwdef struct AkessonModel{MP} <: AbstractModel
    mp::MP
end

function init_akesson_model(mp::ModelParams; kwargs...)

    pars, u0 = init_params(;kwargs...)
    pars = NamedTuple([pair for pair in pars])
    mp = ParametricModels.remake(mp, u0 = u0, p = pars)
    AkessonModel(mp)
end

function (em::AkessonModel)(du, u, p, t)
    eqs!(du, u, p, t)
end