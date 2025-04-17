using Test
using UniqueKronecker: khatri_rao, ⊙

@testset "khatri_rao basic" begin
    A = [1 2; 3 4]
    B = [0 5; 6 7]
    R = khatri_rao(A, B)
    @test size(R) == (4, 2)
    expected = hcat(kron(A[:,1], B[:,1]), kron(A[:,2], B[:,2]))
    @test R == expected

    # Unicode alias
    @test A ⊙ B == R
end

@testset "unary khatri_rao and ⊙" begin
    A = rand(3, 3)
    @test khatri_rao(A) ≈ khatri_rao(A, A)
    @test ⊙(A) ≈ khatri_rao(A, A)
end

@testset "khatri_rao with three matrices" begin
    A = [1 2; 3 4]
    B = [0 5; 6 7]
    C = [1 0; 0 1]
    R3 = khatri_rao(A, B, C)
    expected3 = hcat([kron(kron(A[:,j], B[:,j]), C[:,j]) for j in 1:2]...)
    @test R3 == expected3

    # Alias chaining
    @test A ⊙ B ⊙ C == R3
end

@testset "khatri_rao with integer power" begin
    A = [1 2; 3 4]
    # power 2
    R2 = khatri_rao(A, 2)
    expected2 = hcat([kron(A[:,j], A[:,j]) for j in 1:2]...)
    @test R2 == expected2
    # power 3
    R3 = khatri_rao(A, 3)
    expected3 = hcat([kron(kron(A[:,j], A[:,j]), A[:,j]) for j in 1:2]...)
    @test R3 == expected3

    # Unicode alias for power
    @test A ⊙ 3 == R3
end

@testset "error on mismatched columns" begin
    A = rand(2, 3)
    B = rand(2, 4)
    @test_throws ArgumentError khatri_rao(A, B)
    @test_throws ArgumentError khatri_rao(A, B, A)
end
