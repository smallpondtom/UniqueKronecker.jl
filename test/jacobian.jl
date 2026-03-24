@testset "Jacobian of Unique Kronecker Product" begin

    @testset "Example 3.7 from paper: diffmat_blocks(2, 2)" begin
        Nblocks = diffmat_blocks(2, 2)
        @test length(Nblocks) == 2

        # N_1^{(2)} should be [2 0; 0 1; 0 0]
        N1 = Matrix(Nblocks[1])
        @test N1 == [2 0; 0 1; 0 0]

        # N_2^{(2)} should be [0 0; 1 0; 0 2]
        N2 = Matrix(Nblocks[2])
        @test N2 == [0 0; 1 0; 0 2]
    end

    @testset "diffmat assembles correctly" begin
        N = diffmat(2, 2)
        @test size(N) == (3, 4)   # q2=3, n*q1 = 2*2 = 4
        @test Matrix(N) == [2 0 0 0;
                            0 1 1 0;
                            0 0 0 2]
    end

    @testset "Degree-2 Jacobian matches paper's worked example" begin
        x = [3.0, 5.0]
        J = unique_kronecker_jacobian(x, 2)

        # p_2^u(x) = [x1^2, x1*x2, x2^2] = [9, 15, 25]
        # Jacobian should be:
        #   [2*x1   0  ]     [6  0 ]
        #   [ x2   x1  ]  =  [5  3 ]
        #   [  0   2*x2]     [0  10]
        @test J ≈ [6.0 0.0; 5.0 3.0; 0.0 10.0]
    end

    @testset "Degree-2 specialized matches general" begin
        for n in 2:6
            x = randn(n)
            J_general = unique_kronecker_jacobian(x, 2)
            J_spec = UniqueKronecker.unique_kronecker_jacobian2(x)
            @test J_general ≈ J_spec atol=1e-12
        end
    end

    @testset "Degree-3 specialized matches general" begin
        for n in 2:5
            x = randn(n)
            J_general = unique_kronecker_jacobian(x, 3)
            J_spec = UniqueKronecker.unique_kronecker_jacobian3(x)
            @test J_general ≈ J_spec atol=1e-12
        end
    end

    @testset "Jacobian vs finite differences for degrees 1–5" begin
        ε = 1e-7
        for n in 2:4, deg in 1:5
            x = randn(n)
            J = unique_kronecker_jacobian(x, deg)

            qi = size(J, 1)
            J_fd = zeros(qi, n)
            pu(v) = UniqueKronecker._unique_monomial_vec(v, deg)
            f0 = pu(x)
            for d in 1:n
                xp = copy(x); xp[d] += ε
                xm = copy(x); xm[d] -= ε
                J_fd[:, d] = (pu(xp) - pu(xm)) / (2ε)
            end
            @test J ≈ J_fd atol=1e-5
        end
    end

    @testset "DiffMatCache matches on-the-fly computation" begin
        n = 3
        cache = DiffMatCache(n, 5)
        @test cache.n == 3
        @test cache.pmax == 5

        for deg in 1:5
            x = randn(n)
            J1 = unique_kronecker_jacobian(x, deg)
            J2 = unique_kronecker_jacobian(x, deg, cache)
            @test J1 ≈ J2 atol=1e-14
        end
    end

    @testset "In-place Jacobian" begin
        n = 4; deg = 3
        cache = DiffMatCache(n, deg)
        x = randn(n)
        qi = binomial(n + deg - 1, deg)
        J = zeros(qi, n)
        unique_kronecker_jacobian!(J, x, deg, cache)
        J_ref = unique_kronecker_jacobian(x, deg, cache)
        @test J ≈ J_ref atol=1e-14
    end

    @testset "Matrix (batch) Jacobian" begin
        n = 3; k = 10; deg = 3
        X = randn(n, k)
        cache = DiffMatCache(n, deg)

        Js = unique_kronecker_jacobian(X, deg, cache)
        @test length(Js) == k
        for ℓ in 1:k
            J_ref = unique_kronecker_jacobian(X[:, ℓ], deg, cache)
            @test Js[ℓ] ≈ J_ref atol=1e-14
        end
    end

    @testset "Coupling matrix" begin
        n = 3; p = 3; j = 4
        cache = DiffMatCache(n, p)
        x = randn(n)

        # Create random G matrices
        G = [randn(j, binomial(n + i - 1, i)) for i in 2:p]

        C = UniqueKronecker.coupling_matrix(x, G, cache; pstart=2)
        @test size(C) == (j, n)

        # Verify: C = sum G[i] * J_i
        C_ref = zeros(j, n)
        for (idx, i) in enumerate(2:p)
            Ji = unique_kronecker_jacobian(x, i, cache)
            C_ref += G[idx] * Ji
        end
        @test C ≈ C_ref atol=1e-12
    end

    @testset "Coupling matrix batch" begin
        n = 3; p = 3; j = 4; k = 5
        cache = DiffMatCache(n, p)
        X = randn(n, k)
        G = [randn(j, binomial(n + i - 1, i)) for i in 2:p]

        Cs = UniqueKronecker.coupling_matrix(X, G, cache; pstart=2)
        @test length(Cs) == k
        for ℓ in 1:k
            C_ref = UniqueKronecker.coupling_matrix(X[:, ℓ], G, cache; pstart=2)
            @test Cs[ℓ] ≈ C_ref atol=1e-12
        end
    end

    @testset "Relation to Elimination matrix: L * S * ∇(s⊗i) == ∇(p_i^u)" begin
        for n in 2:4, deg in 2:4
            x = randn(n)
            # Full Kronecker Jacobian via duplication: ∇p_i^u = L * S * ∇(x⊗i)
            # We verify indirectly: D * p_i^u(x) == x⊗i
            D = dupmat(n, deg)
            L = elimat(n, deg)
            S = symmtzrmat(n, deg)
            pu = UniqueKronecker._unique_monomial_vec(x, deg)
            full_kron = reduce(kron, fill(x, deg))
            @test D * pu ≈ full_kron atol=1e-12

            # And verify: L * S * ∇(x⊗i) * e_d matches our Jacobian column d
            # ∇(x⊗i) via formula (44): sum_{α=0}^{i-1} x⊗α ⊗ I ⊗ x⊗(i-1-α)
            Ir = Matrix{Float64}(I, n, n)
            full_jac = zeros(n^deg, n)
            for α in 0:(deg-1)
                left = α == 0 ? ones(1) : reduce(kron, fill(x, α))
                right = (deg - 1 - α) == 0 ? ones(1) : reduce(kron, fill(x, deg - 1 - α))
                full_jac .+= kron(left, kron(Ir, right))
            end
            J_via_dup = Matrix(L) * Matrix(S) * full_jac

            J_ours = unique_kronecker_jacobian(x, deg)
            @test J_via_dup ≈ J_ours atol=1e-10
        end
    end

    @testset "Degree-1 is identity" begin
        for n in 1:5
            x = randn(n)
            J = unique_kronecker_jacobian(x, 1)
            @test J ≈ Matrix{Float64}(I, n, n) atol=1e-14
        end
    end

    @testset "Scalar case (n=1)" begin
        for deg in 1:6
            x = [2.5]
            J = unique_kronecker_jacobian(x, deg)
            @test size(J) == (1, 1)
            # d/dx (x^deg) = deg * x^(deg-1)
            @test J[1, 1] ≈ deg * 2.5^(deg - 1) atol=1e-10
        end
    end

    @testset "Sparsity of differentiation matrices" begin
        for n in [2, 3, 4, 5], i in [2, 3, 4]
            Nblocks = diffmat_blocks(n, i)
            qi = binomial(n + i - 1, i)
            for d in 1:n
                Nd = Nblocks[d]
                # Each row of N_d^{(i)} has at most 1 nonzero entry
                for row in 1:qi
                    @test count(!iszero, Nd[row, :]) ≤ 1
                end
            end
            # Total nnz across all blocks = n * C(n+i-2, i-1)
            total_nnz = sum(nnz(Nblocks[d]) for d in 1:n)
            @test total_nnz == n * binomial(n + i - 2, i - 1)
            # Sum of all values across all blocks = i * q_i
            total_val = sum(sum(Nblocks[d]) for d in 1:n)
            @test total_val == i * qi
        end
    end

end
