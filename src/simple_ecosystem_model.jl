"""
    $SIGNATURES

General multitrophic ecosystem model.

# Args
- `mp`: the model parameters.
- `intinsic_growth_rate`:  the intrinsic growth rate of the population.
- `carrying_capacity`: the maximum population size that can be supported by the environment.
- `competition`:  the effect of population density on the growth rate.
- `resource_conversion_efficiency`: the efficiency of converting resources into population growth.
- `feeding`: feeding rates of the population.

"""
Base.@kwdef struct SimpleEcosystemModel{MP,IGP,CC,CT,RCE,FE} <: AbstractModel
    mp::MP
    intinsic_growth_rate::IGP = default_intinsic_growth_rate
    carrying_capacity::CC = default_carrying_capacity
    competition::CT = default_competition
    resource_conversion_efficiency::RCE = default_resource_conversion_efficiency
    feeding::FE = default_feeding
end

function (model::SimpleEcosystemModel)(du, u, p, t)
    @unpack intinsic_growth_rate, carrying_capacity, competition, feeding, resource_conversion_efficiency = model
    ũ = max.(u, 0.)

    r = intinsic_growth_rate(p, t)
    K = carrying_capacity(p, t)
    A = competition(ũ, p, t)
    ϵ = resource_conversion_efficiency(p, t)
    F = feeding(ũ, p, t)
    
    feed_gains = F * ũ
    pred_loss = F' * ũ
    du .= ũ .*(r .- A ./ K .+ ϵ .* feed_gains .- pred_loss)
end

default_intinsic_growth_rate(p, t) = p.r
default_carrying_capacity(p, t) = p.K
default_competition(u, p, t) = p.A * u
default_resource_conversion_efficiency(p, t) = p.ϵ
default_feeding(u, p, t) = p.q .* p.W ./ (one(eltype(u)) .+ p.q .* p.H .* (p.W * u))