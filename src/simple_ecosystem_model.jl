"""
    $SIGNATURES
    
Based on [Watson et al 2015](https://www.sciencedirect.com/science/article/pii/S0079661114001463) for fish dynamics
and [Steele and Henderson 1981](https://www.jstor.org/stable/2460753) for the nutrient phytoplankton zooplankton dymamics.
"""
Base.@kwdef struct SimpleEcosystemModel{MP,IGP,CC,CT,RCE,FE,HT,AR} <: AbstractModel
    mp::MP
    intinsic_growth_rate::IGP = default_intinsic_growth_rate
    carrying_capacity::CC = default_carrying_capacity
    competition::CT = default_competition
    resource_conversion_efficiency::RCE = default_resource_conversion_efficiency
    feeding::FE = default_feeding
end

function (model::SimpleEcosystemModel)(du, u, p, t)
    @unpack intinsic_growth_rate, carrying_capacity, competition, feeding, handling_time, attack_rate = model
    ũ = max.(u, 0.)

    r = intinsic_growth_rate(p, t)
    K = carrying_capacity(p, t)
    A = competition(ũ, p, t)
    ϵ = resource_conversion_efficiency(p, t)
    F = feeding(ũ, p, t)
    
    # TODO: this is work in progress and should be tested
    # TODO: one should make sure that the sum used for the feeding rates make sense
    du = ũ .*(r .- A ./ K) + ϵ .* sum(F .* N, dims=2) + sum( F .* N, dims=1)
end

default_intinsic_growth_rate(p, t) = p.r
default_carrying_capacity(p, t) = p.K
default_competition(u, p, t) = p.A * u
default_resource_conversion_efficiency(p, t) = p.ϵ
default_feeding(u, p, t) = p.q .* p.W ./ (one(eltype(u)) + p.q .* p.H .* sum(p.W .* u, dims=2))