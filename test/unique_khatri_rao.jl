using Test
using UniqueKronecker: unique_khatri_rao, ⨸, unique_kronecker, unique_kronecker_power

@testset "unique_khatri_rao basic" begin
    A = [1 2; 3 4]
    B = [5 6; 7 8]
    col1 = unique_kronecker(A[:,1], B[:,1])  # [5,7,21]
    col2 = unique_kronecker(A[:,2], B[:,2])  # [12,16,32]
    expected = hcat(col1, col2)
    @test unique_khatri_rao(A, B) == expected
    @test A ⨸ B == expected
    # unary = self‐Khatri–Rao
    @test unique_khatri_rao(A) == unique_khatri_rao(A, A)
    @test ⨸(A) == unique_khatri_rao(A, A)
end

@testset "unique_khatri_rao multiple mats" begin
    A = [1 2; 3 4]
    B = [5 6; 7 8]
    C = [2 0; 0 2]
    cols = [ unique_kronecker(A[:,j], B[:,j], C[:,j]) for j in 1:2 ]
    expected3 = hcat(cols...)
    # This is supposed to be wrong since the unique Kronecker no longer is
    # unique for the third matrix since (A ⊘ A) and A are not symmetric.
    @test unique_khatri_rao(A, B, C) !== expected3
    @test ⨸(A, B, C) !== expected3
end

@testset "unique_khatri_rao power" begin
    A = [1 2; 3 4]
    # d = 2
    cols2 = [ ⊘(A[:,j], 2) for j in 1:2 ]
    expected2 = hcat(cols2...)
    @test unique_khatri_rao(A, 2) == expected2
    @test ⨸(A, 2) == expected2
    # d = 3
    cols3 = [ ⊘(A[:,j], 3) for j in 1:2 ]
    expected3 = hcat(cols3...)
    @test unique_khatri_rao(A, 3) == expected3
end

@testset "unique_khatri_rao errors" begin
    A = rand(2,2)
    C = rand(2,3)
    @test_throws ArgumentError unique_khatri_rao(A, C)
    @test_throws ArgumentError unique_khatri_rao(A, C, A)
    @test_throws ArgumentError unique_khatri_rao(A, 0)
end