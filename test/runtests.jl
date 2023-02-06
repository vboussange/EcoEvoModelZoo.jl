using EcoEvoModelZoo
using Test

@testset "EcoEvoModelZoo.jl" begin
    # Write your tests here.
    include("Huisman_et_al.jl")
    include("McCann_et_al.jl")
    include("Akesson_et_al/model.jl")
    include("Akesson_et_al/init_params.jl")
    include("Akesson_et_al/utils.jl")
end
