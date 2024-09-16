module UniqueKronecker

using LinearAlgebra
using SparseArrays
using Kronecker: ⊗
using StatsBase: countmap
using Combinatorics: permutations, factorial, binomial, with_replacement_combinations

include("unique_kronecker.jl")
include("polytools.jl")
include("vectorize.jl")
include("legacy.jl")

end # module UniqueKronecker
