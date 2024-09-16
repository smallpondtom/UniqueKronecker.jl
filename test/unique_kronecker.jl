@testset "Unique Kronecker" begin
    n = 2
    x = rand(n)
    y = x ⊗ x

    z = x ⊘ x 
    @test all(dupmat(n, 2) * z .== y)

    w = ⊘(x, 3)
    y = x ⊗ x ⊗ x
    @test all(dupmat(n, 3) * w ≈ y)
end
