@testset "Duplication Matrix" begin
    n = 2
    D = UniqueKronecker.dupmat(n, 2)
    @test all(D .== [1 0 0; 0 1 0; 0 1 0; 0 0 1])
end

@testset "Elimination Matrix" begin
    n = 2
    L = UniqueKronecker.elimat(n, 2)
    @test all(L .== [1 0 0 0; 0 1 0 0; 0 0 0 1])
end

@testset "Symmetric Elimination Matrix" begin
    n = 2
    L = UniqueKronecker.elimat(n, 2) * UniqueKronecker.symmtzrmat(n, 2)
    @test all(L .== [1 0 0 0; 0 0.5 0.5 0; 0 0 0 1])
end

@testset "Commutation matrix" begin
    n = 2
    K = UniqueKronecker.commat(n, 2)
    @test all(K .== [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1])
end



