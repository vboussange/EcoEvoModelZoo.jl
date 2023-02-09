using Test
using EcoEvoModelZoo

@testset "init_params" begin
    L = 1; SR = 10; SC = 0
    pars, ic =  EcoEvoModelZoo.init_params_akesson_model(L, SR, SC)
    @test pars isa Dict
    @test ic isa Array
    pars, ic =  EcoEvoModelZoo.init_params_akesson_model(L, SR, SC, width_growth = WidthGrowth{:Standard}())
    @test pars isa Dict
    @test ic isa Array
    pars, ic =  EcoEvoModelZoo.init_params_akesson_model(L, SR, SC, competition = Competition{:Standard}())
    @test pars isa Dict
    @test ic isa Array
    pars, ic =  EcoEvoModelZoo.init_params_akesson_model(L, SR, SR; trophic = Trophic{true}())
    @test pars isa Dict
    @test ic isa Array
end