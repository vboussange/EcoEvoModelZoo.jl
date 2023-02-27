using EcoEvoModelZoo
using ParametricModels
using Random; Random.seed!(1234)
using LinearAlgebra
using UnPack
using Test
using OrdinaryDiffEq
using Statistics
false ? (using PyPlot) : nothing
# Set a random seed for reproduceable behaviour

@testset "Huisman 5 species" begin
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

    p_true = (r = r, m = [m], D = [D], K = K, C = C)

    N0 = [0.1 + i / 100 for i in 1:5]
    R0 = S
    u0 = [N0;R0]
    tsteps = tspan[1]:step:tspan[2]

    model = ResourceCompetition(ModelParams(;p = p_true,
                                            tspan,
                                            u0,
                                            alg,
                                            reltol,
                                            abstol,
                                            saveat = tsteps
                                            ), nN, nR, S)

    # Species 1:5
    sol = simulate(model, u0 = u0)
    @test SciMLBase.successful_retcode(sol.retcode)
    if false # visual test
        _cmap = PyPlot.cm.get_cmap("tab20", 17);
        color_palette = [_cmap(i) for i in 1:17]
        using PyPlot
        fig, ax = plt.subplots()
        ax.plot(tspan[1]:step:tspan[2],Array(sol)')
        # ax.set_yscale("log")
        display(fig)
    end
end

@testset "Smoothmin" begin
    nN = 6
    nR = 3
    r = ones(nN)
    S = [6., 10., 14.]
    K = [1.00  0.90  0.30  1.04  0.34  0.77;
        0.30  1.00 0.90  0.71  1.02  0.76;
        0.90  0.30  1.00  0.46 0.34  1.07]'
    C = [0.04 0.07  0.04  0.10  0.03  0.02
        0.08  0.08  0.10  0.10  0.05  0.17;
        0.14 0.10  0.10  0.16  0.06  0.14;]'
    R = S

    mumin = μ(R, r, K)
    mumean = mean(r .* R' ./ (K .+ R'), dims=2)

    @test all(isapprox.(μsmooth(R, r, K, 0.), mumean, rtol = 0.001))
    @test all(isapprox.(μsmooth(R, r, K, -500.), mumin, rtol = 0.001))
end

@testset "Huisman 5 species with smooth mean" begin
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
    @test SciMLBase.successful_retcode(sol.retcode)
    if false # visual test
        _cmap = PyPlot.cm.get_cmap("tab20", 17);
        color_palette = [_cmap(i) for i in 1:17]
        fig, ax = plt.subplots()
        ax.plot(tspan[1]:step:tspan[2],Array(sol)')
        # ax.set_yscale("log")
        display(fig)
    end
end
