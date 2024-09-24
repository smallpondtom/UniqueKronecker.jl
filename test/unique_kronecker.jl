@testset "Unique Kronecker" begin
    n = 2
    x = rand(n)
    y = x ⊗ x

    z = x ⊘ x 
    @test all(dupmat(n, 2) * z .== y)

    w = ⊘(x, 3)
    y = x ⊗ x ⊗ x
    @test all(dupmat(n, 3) * w ≈ y)

    x = rand(2)
    y = unique_kronecker(x,x)
    z = unique_kronecker(x,x,x)
    w = unique_kronecker(x,x,x,x)
    @test all(⊘(x, x) ≈ y)
    @test all(⊘(x, x, x) ≈ z)
    @test all(⊘(x, x, x, x) ≈ w)

    x = rand(3)
    y = UniqueKronecker.unique_kronecker_power(x,2)
    @test all(⊘(x, 2) ≈ y)

    it = UniqueKronecker.UniqueCombinationIterator(3, 2)
    comb = [1, 1]
    result, _ = iterate(it, comb)
    @test result == [1, 2]
end
