"""
    $SIGNATURES

General multitrophic ecosystem model that models the population dynamics of
multiple species in an ecosystem.

# Arguments
- `mp`: Model parameters of type ModelParams.
- `intinsic_growth_rate`: A function that takes in the model parameters `p` and
  time `t` and returns the intrinsic growth rate of the population.
- `carrying_capacity`: A function that takes in the model parameters `p` and
  time `t` and returns the carrying capacity of the population.
- `competition`: A function that takes in the state `u`, model parameters `p`,
  and time `t` and returns the effect of population density on the growth rate.
- `resource_conversion_efficiency`: A function that takes in the model
  parameters `p` and time `t` and returns the efficiency of converting resources
  into population growth.
- `feeding`: A function that takes in the state `u`, model parameters `p`, and
  time `t` and returns the feeding rates of the population.

# Example 1
Simulating a chaotic 5 compartment food web model proposed in [Post et al.
2000](https://www.jstor.org/stable/177129).

```julia
alg = BS3()
abstol = 1e-6
reltol = 1e-6
tspan = (0., 600)
tsteps = (tspan[1]):1.:(tspan[end])
N = 5 # number of compartment

x_c = 0.15, 
x_p = 0.08
y_c = 2.3
y_p = 1.7
R_0 = 0.25
C_0 = 0.5
ω = 0.2

# constructing simple chain foodweb
# Population vector is organised in such a way
# N = [R1, C1, P, R2, C2]
foodweb = SimpleWeightedDiGraph(N)
add_edge!(foodweb, 2, 1, 1.) # C1 to R1
add_edge!(foodweb, 5, 4, 1.) # C2 to R2
add_edge!(foodweb, 3, 2, ω) # P to C1
add_edge!(foodweb, 3, 5, 1-ω) # P to C2

H = sparse(zeros(N,N))
H[2,1] = 1 / (x_c * y_c); H[5,4] = 1 / (x_c * y_c); 
H[3,2] = 1 / (x_p * y_p); H[3,5] = 1 / (x_p * y_p)

q = sparse(zeros(N,N))
q[2,1] = x_c * y_c / R_0; q[5,4] = x_c * y_c / R_0;
q[3,2] = x_p * y_p / C_0; q[3,5] = x_p * y_p / C_0

p = (r = vcat(1., -x_c, -x_p, 1., -x_c,), 
        K = ones(N), 
        A = diagm(vcat(1,0,0,1,0)), 
        ϵ = ones(N), 
        q = q,
        H = H, 
        W = adjacency_matrix(foodweb))

u0_true = rand(N)


model = SimpleEcosystemModel(;mp = ModelParams(;p,
                                        tspan,
                                        u0 = u0_true,
                                        alg,
                                        reltol,
                                        abstol,
                                        saveat = tsteps,
                                        verbose = false, # suppresses warnings for maxiters
                                        maxiters = 50_000,
                                        ))
sol = simulate(model)
```

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
    T = eltype(u)
    ũ = max.(u, zero(T))

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