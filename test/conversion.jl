@testset "Quadratic matrix conversion" begin    # Settings for the KS equation
    N = 4
    H = zeros(N, N^2)
    for i in 1:N 
        x = rand(N)
        H[i, :] = x ⊗ x
    end
    F = UniqueKronecker.eliminate(H, 2)

    @test all(F .== UniqueKronecker.eliminate(UniqueKronecker.duplicate(F, 2), 2))
    @test all(F .== UniqueKronecker.eliminate(UniqueKronecker.duplicate_symmetric(F, 2), 2))

    H = UniqueKronecker.duplicate(F,2)
    Q = UniqueKronecker.H2Q(H)
    Hnew = Matrix(UniqueKronecker.Q2H(Q))
    @test all(H .== Hnew)
end


@testset "3rd order matrix conversions" begin
    for n in [2, 3, 4, 5, 6, 10]
        x = 1:n # Julia creates Double by default
        x3e = zeros(Int, div(n*(n+1)*(n+2), 6))
        l = 1
        for i = 1:n
            for j = i:n
                for k = j:n
                    x3e[l] = x[i] * x[j] * x[k]
                    l += 1
                end
            end
        end

        L = UniqueKronecker.elimat(n, 3)
        x3 = (x ⊗ x ⊗ x)
        x3e_2 = L * x3
        @test all(x3e_2 .== x3e)

        D = UniqueKronecker.dupmat(n, 3)
        x3_2 = D * x3e
        @test all(x3_2 .== x3)

        G = zeros(n, n^3)
        for i = 1:n
            y = rand(n)
            G[i, :] = y ⊗ y ⊗ y
        end
        E = UniqueKronecker.eliminate(G, 3)
        @test E * x3e ≈ G * x3

        Gs = UniqueKronecker.duplicate_symmetric(E, 3)
        @test Gs * x3 ≈ G * x3

        G2 = UniqueKronecker.duplicate(E, 3)
        @test G2 * x3_2 ≈ G * x3
    end
end

