using Test
using UniqueKronecker: face_split, ⊖

@testset "face_split basic" begin
    A = [1 2; 3 4]
    B = [5 6; 7 8]
    # row1: kron([1,2],[5,6]) = [5,6,10,12]
    # row2: kron([3,4],[7,8]) = [21,24,28,32]
    expected = [5 6 10 12;
                21 24 28 32]
    @test face_split(A, B) == expected
    @test A ⊖ B == expected
end

@testset "face_split unary" begin
    A = rand(3,4)
    @test face_split(A) ≈ face_split(A, A)
    @test ⊖(A) ≈ face_split(A, A)
end

@testset "face_split multiple mats" begin
    A = [1 2; 3 4]
    B = [5 6; 7 8]
    C = [2 3; 4 5]
    # row1: kron(kron([1,2],[5,6]),[2,3])
    row1 = kron(kron([1,2], [5,6]), [2,3])
    # row2: kron(kron([3,4],[7,8]),[4,5])
    row2 = kron(kron([3,4], [7,8]), [4,5])
    expected3 = vcat(row1', row2')
    @test face_split(A, B, C) == expected3
    @test ⊖(A, B, C) == expected3
end

@testset "face_split power" begin
    A = [1 2; 3 4]
    # d = 2
    r1_2 = kron([1,2], [1,2])
    r2_2 = kron([3,4], [3,4])
    expected2 = vcat(r1_2', r2_2')
    @test face_split(A, 2) == expected2
    @test A ⊖ 2 == expected2

    # d = 3
    r1_3 = reduce(kron, ntuple(_->[1,2], 3))
    r2_3 = reduce(kron, ntuple(_->[3,4], 3))
    expected3 = vcat(r1_3', r2_3')
    @test face_split(A, 3) == expected3
    @test A ⊖ 3 == expected3
end

@testset "face_split errors" begin
    A = rand(2,3)
    B = rand(3,3)
    @test_throws ArgumentError face_split(A, B)
    @test_throws ArgumentError face_split(A, 0)
end