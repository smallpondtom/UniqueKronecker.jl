@testset "vech" begin
    A = [1 2; 3 4]
    a = [1; 3; 4]
    @test all(a .== UniqueKronecker.vech(A))
end


@testset "inverse vectorization" begin
    A = [1 2; 3 4]
    a = vec(A)
    @test all(A .== UniqueKronecker.invec(a, 2, 2))
end
