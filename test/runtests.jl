using LinearAlgebra
using SparseArrays
using UniqueKronecker
using Test
using Kronecker 
using Random

function testfile(file, testname=defaultname(file))
    println("running test file $(file)")
    @testset "$testname" begin; include(file); end
    return
end
defaultname(file) = uppercasefirst(replace(splitext(basename(file))[1], '_' => ' '))

@testset "UniqueKronecker" begin
    testfile("unique_kronecker.jl")
    testfile("circulant_kronecker.jl")
    testfile("special_matrices.jl")
    testfile("khatri_rao.jl")
    testfile("unique_khatri_rao.jl")
    testfile("face_split.jl")
    testfile("unique_face_split.jl")
    testfile("create_poly.jl")
    testfile("conversion.jl")
    testfile("snapshot.jl")
    testfile("vectorize.jl")
    testfile("jacobian.jl")
    testfile("legacy.jl")
end