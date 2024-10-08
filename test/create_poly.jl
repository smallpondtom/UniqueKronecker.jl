# @testset "Making polynomial operator" begin
#     idx = [(1, 1, 1), (1, 1, 2)]
#     val = [1.0, 2.0]
#     H = UniqueKronecker.makeQuadOp(3, idx, val, which_quad_term="H", symmetric=true)
#     F = UniqueKronecker.makeQuadOp(3, idx, val, which_quad_term="F", symmetric=false)

#     A2 = UniqueKronecker.make_poly_op(3, idx, val, nonredundant=false, symmetric=true)
#     A2u = UniqueKronecker.make_poly_op(3, idx, val, nonredundant=true, symmetric=false)
#     @test all(H .== A2)
#     @test all(F .== A2u)
# end

@testset "Making polynomial operator" begin
    N = 3
    K = 3
    for n in 2:N
        for k in 2:K
            du = binomial(n+k-1, k)
            idx = Array[]
            val = Float64[]

            for i in 1:du
                push!(idx, rand(1:n, k+1))
                push!(val, rand(-2:2))
            end
            idx = map(x -> Tuple(x), idx)

            if k == 2
                H = UniqueKronecker.makeQuadOp(n, idx, val, which_quad_term="H", symmetric=true)
                F = UniqueKronecker.makeQuadOp(n, idx, val, which_quad_term="F", symmetric=false)
            elseif k == 3
                H = UniqueKronecker.makeCubicOp(n, idx, val, which_cubic_term="G", symmetric=true)
                F = UniqueKronecker.makeCubicOp(n, idx, val, which_cubic_term="E", symmetric=false)
            end

            Ak_1 = UniqueKronecker.make_poly_op(n, idx, val, nonredundant=false, symmetric=true)
            Aku_1 = UniqueKronecker.make_poly_op(n, idx, val, nonredundant=true, symmetric=false)
            if k == 2 || k == 3
                @test all(H .== Ak_1)
                @test all(F .== Aku_1)
            end

            # Ak_2 = UniqueKronecker.make_poly_op_faster(n, idx, val, nonredundant=false, symmetric=true)
            # Aku_2 = UniqueKronecker.make_poly_op_faster(n, idx, val, nonredundant=true, symmetric=false)
            # @test all(Ak_1 .== Ak_2)
            # @test all(Aku_1 .== Aku_2)

            # Ak = UniqueKronecker.make_poly_op_parallel(n, idx, val, nonredundant=false, symmetric=true)
            # Aku = UniqueKronecker.make_poly_op_parallel(n, idx, val, nonredundant=true, symmetric=false)
            # @test all(Ak_1 .== Ak)
            # @test all(Aku_1 .== Aku)
        end
    end
end
