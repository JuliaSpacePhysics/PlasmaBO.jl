using PlasmaBO
using Test
using Aqua

@testset "PlasmaBO.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(PlasmaBO)
    end
    # Write your tests here.
end
