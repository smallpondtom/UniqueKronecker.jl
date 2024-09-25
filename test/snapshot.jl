@testset "Creating polynomial snapshot matrices" begin 
    n = 3
    K = 10
    X = rand(n, K)
    X2 = kron_snapshot_matrix(X, 2)
    X2u = unique_kron_snapshot_matrix(X, 2)

    @test all(elimat(n,2) * X2 .== X2u)
end

@testset "Circulant Kronecker snapshot matrix" begin
    X1 = [1 2; 3 4]
    X2 = [5 6; 7 8]
    X3 = [9 10; 11 12]

    X = circulant_kron_snapshot_matrix(X1, X2, X3)
    X_nr = circulant_kron_snapshot_matrix(X1, X2, X3; redundant=false)
    L23 = elimat(2, 3)

    for i in 1:2
        x1 = X1[:,i]
        x2 = X2[:,i]
        x3 = X3[:,i]
        x = circulant_kronecker(x1, x2, x3)
        x_nr = L23 * vec(x)
        @test X[:,i] == vec(x)
        @test X_nr[:,i] == vec(x_nr)
    end
end
