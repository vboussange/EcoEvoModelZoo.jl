var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = EcoEvoModelZoo","category":"page"},{"location":"#EcoEvoModelZoo","page":"Home","title":"EcoEvoModelZoo","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"DocStringExtensions.README","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [EcoEvoModelZoo]\nPrivate = false\n<!-- Filter = t -> typeof(t) <: AbstractModel -->","category":"page"},{"location":"#EcoEvoModelZoo.AkessonModel","page":"Home","title":"EcoEvoModelZoo.AkessonModel","text":"This model is inspired from Akesson et al. 2021.\n\nSpecialized variants are provided in other AbstractModels.\n\n\n\n\n\n","category":"type"},{"location":"#EcoEvoModelZoo.EcoEvoGraph","page":"Home","title":"EcoEvoModelZoo.EcoEvoGraph","text":"This model is inspired from Boussange & Pellissier. 2022\n\nArguments\n\nmp: the model parameters\ng: the spatial graph\nphen_space: the discretized phenotypic space\nbirth_fn: the birth function. Should be of the form birth_fn(x, s, p)\ncompetition_fn: the competition function. Should be of the form competition_fn(u_xs2, x, s1, s2, p)\n\nMathematical model\n\n∂ₜuˣ(s) = uˣ(s)[birth_fn(x, s, p) - uˣ(s)∫competition_fn(uˣ(s₂), x, s, s₂, p) ds₂] + 1/2 μ σ² Δₛuˣ(s) + m L u(s)\n\nwhere L is the Laplacian matrix of the spatial graph g.\n\nExample 1\n\nusing ParametricModels\nusing UnPack\nusing DocStringExtensions\nusing Statistics\nusing ComponentArrays\nusing Distributions\nusing OrdinaryDiffEq\nusing UnPack\n## Defining phenotypic space and spatial graph\nM = 7\ng = star_graph(M)\n\nrS = 1f0\ndS = 0.02f0\nphen_space = collect(range(-rS,rS,step=dS)) #grid\n\n##\n# Defining birth and competition functions\nsoptim = 0.5f0 * [-1,1,-1,1,-1,1,-1]\nbirth_fn(x, s, p) = max(0f0, 1f0 - (soptim[x] - s)^2)\ncompetition_fn(u_xs2, x, s1, s2, p) = u_xs2 ./ p.K\n\n## Defining parameters\nσ_mu = 5f-2;\nmu = 0.1f0\nm = 0.1\nK = 1.\n\np = ComponentArray(σ_mu = σ_mu,\n                    mu = mu,\n                    m = m,\n                    K = K)\n\n\n## rest of the simulation\ntend = 1000f0\ntspan = (0f0,tend)\ntsteps = (tspan[1]):1.:(tspan[end])\nu0 = vcat([K .* pdf.(Normal(so,σ_mu),phen_space') for so in soptim]...)\n\nmp = ModelParams(;p,\n                tspan,\n                u0,\n                alg=Tsit5(),\n                saveat = tsteps)\n\nmodel = EcoEvoGraph(mp, g, phen_space, birth_fn, competition_fn)\n\n@time sol = simulate(model)\n\n# Plotting results!\nusing PythonCall; plt = pyimport(\"matplotlib.pyplot\")\nfig, ax = plt.subplots(1)\nax.plot(model.phen_space, sol[end][1,:], color = \"tab:red\", label = \"Node 1\")\nax.plot(model.phen_space, sol[end][2,:], color = \"tab:blue\", label = \"Node 2\")\nax.set_xlabel(\"Phenotype\")\nax.set_ylabel(\"Population number\")\ndisplay(fig)\n\nExample 2: trait dependent competition\n\n# Trait dependent compeition\nrS = 3f0\ndS = 0.02f0\nphen_space = collect(range(-rS,rS,step=dS)) #grid\n\np = ComponentArray(σ_mu = σ_mu,\n                    mu = mu,\n                    m = m,\n                    K = K,\n                    σ_α = 0.5)\n# trait dependent competition\ncompetition_fn(u_xs2, x, s1, s2, p) = u_xs2 * exp(- 0.5f0 * (s1 - s2)^2 ./ p.σ_α^2) ./ p.K\n\nu0 = vcat([K .* pdf.(Normal(so,σ_mu),phen_space') for so in soptim]...)\n\nmp = ModelParams(;p,\n                tspan,\n                u0,\n                alg=Tsit5(),\n                saveat = tsteps)\nmodel = EcoEvoGraph(mp, g, phen_space, birth_fn, competition_fn)\n\nsol = simulate(model)\nfig, ax = plt.subplots(1)\nax.set_title(\"Trait-dependent competition\")\nax.plot(model.phen_space, sol[end][1,:], color = \"tab:red\", label = \"Node 1\")\nax.plot(model.phen_space, sol[end][2,:], color = \"tab:blue\", label = \"Node 2\")\nax.set_xlabel(\"Phenotype\")\nax.set_ylabel(\"Population number\")\ndisplay(fig)\n\n\n\n\n\n","category":"type"},{"location":"#EcoEvoModelZoo.EcosystemModelMcCann","page":"Home","title":"EcoEvoModelZoo.EcosystemModelMcCann","text":"EcosystemModelMcKann\n\nThis model is inspired from McCann 1994.  Similar to EcosystemModelOmnivory, but without monivory.\n\nModel parameters\n\nx_c, x_p: mass-specific metabolic rate of consumers and predators\ny_c, y_p: ingestion rate per unit metabolic rate of consumers and predators.\nR_0, C_0: half saturation densities for the type II functional responses of the consumers and predators\n\nExample\n\n```julia alg = BS3() abstol = 1e-6 reltol = 1e-6 tspan = (0., 800) tsteps = 550:4:800 ptrue = (xc = [0.4],          xp = [0.08],          yc = [2.01],         yp = [5.00],          R0 = [0.16129],          C_0 = [0.5])\n\nu0_true = [0.5,0.8,0.5]\n\nmodel = EcosystemModelMcCann(ModelParams(;p = ptrue,                                         tspan,                                         u0 = u0true,                                         alg,                                         reltol,                                         abstol,                                         saveat = tsteps,                                         verbose = false, # suppresses warnings for maxiters                                         maxiters = 50000,                                         )) sol = simulate(model, u0 = u0true)\n\n\n\n\n\n","category":"type"},{"location":"#EcoEvoModelZoo.EcosystemModelOmnivory","page":"Home","title":"EcoEvoModelZoo.EcosystemModelOmnivory","text":"EcosystemModelOmnivory\n\nThis model is inspired from McCann 1997.\n\nModel parameters\n\nx_c, x_p: mass-specific metabolic rate of consumers and predators\ny_c, y_pr, y_pc: ingestion rate per unit metabolic rate of consumers and predators.\nR_0, R_02, C_0 half saturation densities for the type II functional responses of the consumers and predators\nω: omnivory strength\n\nExample\n\nalg = BS3()\nabstol = 1e-6\nreltol = 1e-6\ntspan = (0., 800)\ntsteps = 550:4:800\n\np_true = (x_c = [0.4], \n        x_p = [0.08], \n        y_c = [2.01], \n        y_pr = [2.00], \n        y_pc = [5.0], \n        R_0 = [0.16129], \n        R_02 = [ 0.5], \n        C_0 = [0.5],\n        ω =[ 0.4])\n\nu0_true = [0.5,0.8,0.5]\n\n\nmodel = EcosystemModelOmnivory(ModelParams(;p = p_true,\n                                        tspan,\n                                        u0 = u0_true,\n                                        alg,\n                                        reltol,\n                                        abstol,\n                                        saveat = tsteps,\n                                        verbose = false, # suppresses warnings for maxiters\n                                        maxiters = 50_000,\n                                        ))\nsol = simulate(model, u0 = u0_true)\n\n```\n\n\n\n\n\n","category":"type"},{"location":"#EcoEvoModelZoo.Landscape","page":"Home","title":"EcoEvoModelZoo.Landscape","text":"Returns landscape parameters\n\n\n\n\n\n","category":"type"},{"location":"#EcoEvoModelZoo.ResourceCompetition","page":"Home","title":"EcoEvoModelZoo.ResourceCompetition","text":"This model is inspired from Huisman et al. 1999 Nature.\n\n\n\n\n\n","category":"type"},{"location":"#EcoEvoModelZoo.ResourceCompetitionSmoothMin","page":"Home","title":"EcoEvoModelZoo.ResourceCompetitionSmoothMin","text":"This model is inspired from Huisman et al. 1999 Nature.,  but where Leibig's law is replaced by  imperfect substituable resources (smooth minimum). The smooth min function is parametrized by  s, which is a trainable parameter.\n\n\n\n\n\n","category":"type"},{"location":"#EcoEvoModelZoo.SimpleEcosystemModel","page":"Home","title":"EcoEvoModelZoo.SimpleEcosystemModel","text":"General multitrophic ecosystem model.\n\nArgs\n\nmp: the model parameters.\nintinsic_growth_rate:  the intrinsic growth rate of the population.\ncarrying_capacity: the maximum population size that can be supported by the environment.\ncompetition:  the effect of population density on the growth rate.\nresource_conversion_efficiency: the efficiency of converting resources into population growth.\nfeeding: feeding rates of the population.\n\n\n\n\n\n","category":"type"},{"location":"#EcoEvoModelZoo.Temperature","page":"Home","title":"EcoEvoModelZoo.Temperature","text":"Temperature as a function of space (x), time (t), and some climate parameters\n\n\n\n\n\n","category":"type"},{"location":"#EcoEvoModelZoo.funcresp!-Tuple{Any, Any, Any, Trophic{true}}","page":"Home","title":"EcoEvoModelZoo.funcresp!","text":"funcresp!(F, n, p, _)\n\n\nType II functional response Input:\n\nn: Vector of population densities of all species in a given patch\nTh: Vector of handling times (with dummy values for resource species)\narate: Vector of attack rates (with dummy values for resource species)\nW: Adjacency matrix of trophic network; W(i,j)=1 if i eats j and 0 otherwise\n\nOutput:\n\nA matrix F(i,j), the feeding rate of consumer i on resource j\n\n\n\n\n\n","category":"method"},{"location":"#EcoEvoModelZoo.generate_network-Tuple{Int64, Int64}","page":"Home","title":"EcoEvoModelZoo.generate_network","text":"generate_network(SR::Int, SC::Int)\n\nreturn matrix W[i,j], which is nonzero if consumer i eats resource j SR is the number of resource, SC the number of consumer species.\n\nThe bipartite network is generated as follows.  First, both resources and consumers are labeled consecutively,  based on their initial temperature adaptations:  resource 1 / consumer 1 are the most cold-adapted,  and resource S / consumer S the most warm-adapted.  Next, we always put a feeding link between consumer i and resource i.  Finally, each consumer is randomly linked to four other resource species.  This yields a feeding network where every consumer is  connected to five resources altogether.\n\nNote: it seems that the definition of W_ij is that it determines which resource i is eat by consumer j\n\n\n\n\n\n","category":"method"},{"location":"#EcoEvoModelZoo.smoothstep-Tuple{Any}","page":"Home","title":"EcoEvoModelZoo.smoothstep","text":"Apply twice continuously differentiable smoothed step function to a number x Input: x: Distance from pole, measured in units of the pole-to-equator distance Output: 0 if x < 0; 10x^3-15x^4+6*x^5 if 0 <= x <= 1; otherwise 1\n\n\n\n\n\n","category":"method"}]
}
