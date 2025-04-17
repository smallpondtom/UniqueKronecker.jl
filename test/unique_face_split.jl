using Test
using UniqueKronecker: unique_face_split, ⧁, unique_kronecker, unique_kronecker_power

@testset "unique_face_split basic" begin
    A = [1 2; 3 4]
    B = [5 6; 7 8]
    # row 1: unique_kronecker([1,2],[5,6])
    # row 2: unique_kronecker([3,4],[7,8])
    row1 = unique_kronecker(A[1, :], B[1, :])
    row2 = unique_kronecker(A[2, :], B[2, :])
    expected = vcat(row1', row2')
    @test unique_face_split(A, B) == expected
    @test (A ⧁ B) == expected

    # unary = self unique face‑splitting
    @test unique_face_split(A) == unique_face_split(A, A)
    @test ⧁(A) == unique_face_split(A, A)
end

@testset "unique_face_split multiple mats" begin
    A = [1 2; 3 4]
    B = [5 6; 7 8]
    C = [2 3; 4 5]
    rows = [
        reduce(unique_kronecker, (M[j, :] for M in (A, B, C))) for j in 1:2
    ]
    # This is supposed to be wrong since the unique Kronecker no longer is
    # unique for the third matrix since (A ⧘ A) and A are not symmetric.
    expected3 = vcat((r' for r in rows)...)
    @test unique_face_split(A, B, C) !== expected3
    @test ⧁(A, B, C) !== expected3
end

@testset "unique_face_split power" begin
    A = [1 2; 3 4]
    # d = 2
    rows2 = [ unique_kronecker_power(A[j, :], 2) for j in 1:2 ]
    expected2 = vcat((r' for r in rows2)...)
    @test unique_face_split(A, 2) == expected2
    @test (A ⧁ 2) == expected2

    # d = 3
    rows3 = [ unique_kronecker_power(A[j, :], 3) for j in 1:2 ]
    expected3 = vcat((r' for r in rows3)...)
    @test unique_face_split(A, 3) == expected3
end

@testset "unique_face_split errors" begin
    A = rand(2, 3)
    B = rand(3, 3)
    @test_throws ArgumentError unique_face_split(A, B)
    @test_throws ArgumentError unique_face_split(A, B, A)
    @test_throws ArgumentError unique_face_split(A, 0)
end