using ParametricModels
using UnPack
using DocStringExtensions

"""
    EcosystemModelOmnivory

This model is inspired from [McCann 1997](10.1098/rspb.1997.0172).

# Model parameters
- `x_c, x_p`: mass-specific metabolic rate of consumers and predators
- `y_c, y_pr, y_pc`: ingestion rate per unit metabolic rate of consumers and predators.
-  `R_0, R_02, C_0` half saturation densities for the type II functional responses of the consumers and predators
- `ω`: omnivory strength
# Example
```julia
alg = BS3()
abstol = 1e-6
reltol = 1e-6
tspan = (0., 800)
tsteps = 550:4:800

p_true = (x_c = [0.4], 
        x_p = [0.08], 
        y_c = [2.01], 
        y_pr = [2.00], 
        y_pc = [5.0], 
        R_0 = [0.16129], 
        R_02 = [ 0.5], 
        C_0 = [0.5],
        ω =[ 0.4])

u0_true = [0.5,0.8,0.5]


model = EcosystemModelOmnivory(ModelParams(;p = p_true,
                                        tspan,
                                        u0 = u0_true,
                                        alg,
                                        reltol,
                                        abstol,
                                        saveat = tsteps,
                                        verbose = false, # suppresses warnings for maxiters
                                        maxiters = 50_000,
                                        ))
sol = simulate(model, u0 = u0_true)
```
```
"""
Base.@kwdef struct EcosystemModelOmnivory{MP} <: AbstractModel
    mp::MP
end

function (em::EcosystemModelOmnivory)(du, u, p, t)
    ũ = max.(u, 0.)
    @unpack x_c, x_p, y_c, y_pr, y_pc, R_0, R_02, C_0, ω = p
    R, C, P = ũ
    du[1] = R * (1. - R) - x_c[] * y_c[] * C * R / (R + R_0[]) - ω[] * x_p[] * y_pr[] * P * R / (R_02[] + (1. - ω[]) * C + ω[] * R) #R
    du[2] = - x_c[] * C * ( 1. - y_c[] * R / (R + R_0[])) - (1. - ω[]) * x_p[] * y_pc[] * P * C / (ω[] * R + (1. - ω[]) * C + C_0[]) #C
    du[3] = - x_p[] * P + (1. - ω[]) * x_p[] * y_pc[] * C * P / (ω[] * R + (1. - ω[]) * C + C_0[]) + ω[] * x_p[] * y_pr[] * P * R / (ω[] * R + (1. - ω[]) * C + R_02[]) #P
end


"""
    EcosystemModelMcKann

This model is inspired from [McCann 1994](http://doi.wiley.com/10.2307/1939558). 
Similar to EcosystemModelOmnivory, but without monivory.

# Model parameters
- `x_c, x_p`: mass-specific metabolic rate of consumers and predators
- `y_c, y_p`: ingestion rate per unit metabolic rate of consumers and predators.
-  `R_0, C_0`: half saturation densities for the type II functional responses of the consumers and predators

# Example
```julia
alg = BS3()
abstol = 1e-6
reltol = 1e-6
tspan = (0., 800)
tsteps = 550:4:800
p_true = (x_c = [0.4], 
        x_p = [0.08], 
        y_c = [2.01],
        y_p = [5.00], 
        R_0 = [0.16129], 
        C_0 = [0.5])

u0_true = [0.5,0.8,0.5]


model = EcosystemModelMcCann(ModelParams(;p = p_true,
                                        tspan,
                                        u0 = u0_true,
                                        alg,
                                        reltol,
                                        abstol,
                                        saveat = tsteps,
                                        verbose = false, # suppresses warnings for maxiters
                                        maxiters = 50_000,
                                        ))
sol = simulate(model, u0 = u0_true)
"""
Base.@kwdef struct EcosystemModelMcCann{MP} <: AbstractModel
    mp::MP # parameter vector
end

function (em::EcosystemModelMcCann)(du, u, p, t) 
    ũ = max.(u, 0.)
    @unpack x_c, x_p, y_c, y_p, R_0, C_0 = p
    R, C, P = ũ
    du[1] = R * (1. - R) - x_c[] * y_c[] * C * R / (R + R_0[]) #R
    du[2] = - x_c[] * C * ( 1. - y_c[] * R / (R + R_0[])) - x_p[] * y_p[] * P * C / (C + C_0[]) #C
    du[3] = - x_p[] * P + x_p[] * y_p[] * C * P / (C + C_0[]) #P
end