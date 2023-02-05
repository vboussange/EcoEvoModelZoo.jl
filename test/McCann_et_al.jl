using ParametricModels
using UnPack
using DocStringExtensions

"""
    EcosystemModelOmnivory

This model is inspired from [McCann 1997](10.1098/rspb.1997.0172).
"""
Base.@kwdef struct EcosystemModelOmnivory{MP} <: AbstractModel
    mp::MP
end

function (em::EcosystemModelOmnivory)(du, u, p, t)
    ũ = max.(u, 0.)
    @unpack x_c, x_p, y_c, y_pr, y_pc, R_0, R_02, C_0, ω = p
    R, C, P = ũ
    du[1] = R * (1. - R) - x_c * y_c * C * R / (R + R_0) - ω * x_p * y_pr * P * R / (R_02 + (1. - ω) * C + ω * R) #R
    du[2] = - x_c * C * ( 1. - y_c * R / (R + R_0)) - (1. - ω) * x_p * y_pc * P * C / (ω * R + (1. - ω) * C + C_0) #C
    du[3] = - x_p * P + (1. - ω) * x_p * y_pc * C * P / (ω * R + (1. - ω) * C + C_0) + ω * x_p * y_pr * P * R / (ω * R + (1. - ω) * C + R_02) #P
end


"""
    EcosystemModelMcKann

This model is inspired from [McCann 1994](http://doi.wiley.com/10.2307/1939558). 
Similar to EcosystemModelOmnivory, but without monivory.
"""
Base.@kwdef struct EcosystemModelMcKann{MP} <: AbstractModel
    mp::MP # parameter vector
end

function (em::EcosystemModelMcKann)(du, u, p, t) 
    ũ = max.(u, 0.)
    @unpack x_c, x_p, y_c, y_p, R_0, C_0 = p
    R, C, P = ũ
    du[1] = R * (1. - R) - x_c[] * y_c[] * C * R / (R + R_0[]) #R
    du[2] = - x_c[] * C * ( 1. - y_c[] * R / (R + R_0[])) - x_p[] * y_p[] * P * C / (C + C_0[]) #C
    du[3] = - x_p[] * P + x_p[] * y_p[] * C * P / (C + C_0[]) #P
end