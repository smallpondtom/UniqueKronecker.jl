@testset "circulant Kronecker" begin
    # Example usage
    # Two vectors
    x = [1, 2]
    y = [3, 4]

    result_xy = x ⊛ y
    result_xy2 = circulant_kronecker(x, y)
    @test result_xy ≈ [6, 10, 10, 16]
    @test result_xy2 ≈ [6, 10, 10, 16]

    # Three vectors
    z = [5, 6]
    result_xyz = ⊛(x, y, z)
    result_xyz2 = circulant_kronecker(x, y, z)
    @test result_xyz ≈ [45, 68, 68, 100, 68, 100, 100, 144]
    @test result_xyz2 ≈ [45, 68, 68, 100, 68, 100, 100, 144]

    # Matrices
    A = [1 2; 3 4]
    B = [5 6; 7 8]
    C = [9 10; 11 12]

    result_ABC = ⊛(A, B,  C)
    result_ABC2 = circulant_kronecker(A, B, C)
    @test result_ABC ≈ [
        135  194  194  268  194  268  268   360
        253  312  342  416  342  416  452   544
        253  342  312  416  342  452  416   544
        431  520  520  624  562  672  672   800
        253  342  342  452  312  416  416   544
        431  520  562  672  520  624  672   800
        431  562  520  672  520  672  624   800
        693  824  824  976  824  976  976  1152
    ]
    @test result_ABC2 ≈ [
        135  194  194  268  194  268  268   360
        253  312  342  416  342  416  452   544
        253  342  312  416  342  452  416   544
        431  520  520  624  562  672  672   800
        253  342  342  452  312  416  416   544
        431  520  562  672  520  624  672   800
        431  562  520  672  520  672  624   800
        693  824  824  976  824  976  976  1152
    ]

    x = [1,2]
    y = [3,4]
    z = [5,6]
    w = [7,8]
    result_xyzw = circulant_kronecker(x, y, z, w)
    @test size(result_xyzw) == (16,1)
end