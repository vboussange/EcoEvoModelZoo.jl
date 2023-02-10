using Test
using EcoEvoModelZoo

@testset "init_params" begin
    L = 1; SR = 10; SC = 0
    land = Landscape(L)
    temp = Temperature()

    pars, ic =  EcoEvoModelZoo.init_params_akesson_model(land, temp, SR, SC)
    @test pars isa Dict
    @test ic isa Array
    pars, ic =  EcoEvoModelZoo.init_params_akesson_model(land, temp, SR, SC, width_growth = WidthGrowth{:Standard}())
    @test pars isa Dict
    @test ic isa Array
    pars, ic =  EcoEvoModelZoo.init_params_akesson_model(land, temp, SR, SC, competition = Competition{:Standard}())
    @test pars isa Dict
    @test ic isa Array
    pars, ic =  EcoEvoModelZoo.init_params_akesson_model(land, temp, SR, SR; trophic = Trophic{true}())
    @test pars isa Dict
    @test ic isa Array
end