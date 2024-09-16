@testset "Making polynomial operator" begin
    idx = [(1, 1, 1), (1, 1, 2)]
    val = [1.0, 2.0]
    H = UniqueKronecker.makeQuadOp(3, idx, val, which_quad_term="H", symmetric=true)
    F = UniqueKronecker.makeQuadOp(3, idx, val, which_quad_term="F", symmetric=false)

    A2 = UniqueKronecker.makePolyOp(3, idx, val, nonredundant=false, symmetric=true)
    A2u = UniqueKronecker.makePolyOp(3, idx, val, nonredundant=true, symmetric=false)
    @test all(H .== A2)
    @test all(F .== A2u)
end
