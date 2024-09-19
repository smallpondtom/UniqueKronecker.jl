@testset "Run legacy code 1" begin
    H = [
        0.0  1.0  0.0  0.0  1.0  0.0  0.0  0.0  1.0
        0.0  0.0  1.0  0.0  1.0  0.0  0.0  0.0  1.0
        1.0  0.0  0.0  0.0  1.0  1.0  0.0  0.0  0.0
    ]
    H2 = UniqueKronecker.insert2H(H, 4)
    H2test = [
        0.0  1.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0
        0.0  0.0  1.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0
        1.0  0.0  0.0  0.0  0.0  1.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
        0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
    ]
    @test all(H2 .== H2test)
    @test all(H .== UniqueKronecker.extractH(H2, 3))

    F = UniqueKronecker.eliminate(H, 2)
    F2 = UniqueKronecker.insert2F(F, 4)
    F2test = [
        0.0  1.0  0.0  0.0  1.0  0.0  0.0  1.0  0.0  0.0
        0.0  0.0  1.0  0.0  1.0  0.0  0.0  1.0  0.0  0.0
        1.0  0.0  0.0  0.0  1.0  1.0  0.0  0.0  0.0  0.0
        0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
    ]
    @test all(F2 .== F2test)
    @test all(F .== UniqueKronecker.extractF(F2, 3))
end

@testset "Run legacy code 2" begin
    n, j, k = 4, 2, 3
    X = rand(2,2)
    BL = UniqueKronecker.insert2bilin(X, 3, 1)

    idx = [(1, 1, 1), (1, 1, 2)]
    val = [1.0, 2.0]
    H = UniqueKronecker.makeQuadOp(3, idx, val, which_quad_term="H", symmetric=true)
    F = UniqueKronecker.makeQuadOp(3, idx, val, which_quad_term="F", symmetric=false)
    Q = UniqueKronecker.makeQuadOp(3, idx, val, which_quad_term="Q", symmetric=false)
    @test all(H .== UniqueKronecker.duplicate_symmetric(F, 2))
    @test all(H .== UniqueKronecker.Q2H(Q))

    idx = [(1, 1, 1, 1), (1, 1, 1, 2)]
    val = [1.0, 2.0]
    G = UniqueKronecker.makeCubicOp(3, idx, val, which_cubic_term="G", symmetric=true)
    E = UniqueKronecker.makeCubicOp(3, idx, val, which_cubic_term="E", symmetric=false)
    @test all(G .== UniqueKronecker.duplicate_symmetric(E, 3))

    G = UniqueKronecker.makeIdentityCubicOp(3, which_cubic_term="G")
    E = UniqueKronecker.makeIdentityCubicOp(3, which_cubic_term="E")
    @test size(G) == (3, 27)
    @test size(E) == (3, 10)
end
