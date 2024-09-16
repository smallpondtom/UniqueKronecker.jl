@testset "Creating polynomial snapshot matrices" begin 
    n = 3
    K = 10
    X = rand(n, K)
    X2 = UniqueKronecker.kron_snapshot_matrix(X, 2)
    X2u = UniqueKronecker.unique_kron_snapshot_matrix(X, 2)

    @test all(elimat(n,2) * X2 .== X2u)
end
