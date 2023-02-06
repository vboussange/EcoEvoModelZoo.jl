using Test

@testset "init_params" begin
    EcoEvoModelZoo.init_params() isa Dict
end