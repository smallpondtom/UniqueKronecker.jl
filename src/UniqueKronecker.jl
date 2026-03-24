module UniqueKronecker

using LinearAlgebra
using SparseArrays
using Kronecker: ⊗, kronecker
using StatsBase: countmap
using Combinatorics: permutations, factorial, binomial, with_replacement_combinations

include("unique_kronecker.jl")
include("circulant_kronecker.jl")
include("khatri_rao.jl")
include("unique_khatri_rao.jl")
include("face_split.jl")
include("unique_face_split.jl")
include("polytools.jl")
include("vectorize.jl")
include("jacobian.jl")
include("legacy.jl")

end # module UniqueKronecker
