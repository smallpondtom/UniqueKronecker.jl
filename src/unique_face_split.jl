export unique_face_split, ⧁

"""
    unique_face_split(A::AbstractMatrix, B::AbstractMatrix)

Row-wise "unique face-splitting" product of A and B: for each row j,
compute `unique_kronecker(A[j, :], B[j, :])` and stack the transposed
results as rows of the output.
"""
function unique_face_split(A::AbstractMatrix, B::AbstractMatrix)
    size(A, 1) == size(B, 1) ||
      throw(ArgumentError("matrices must have same number of rows"))
    return vcat(map(transpose ∘ unique_kronecker, eachrow(A), eachrow(B))...)
end

unique_face_split(A::AbstractMatrix) = unique_face_split(A, A)

"""
    unique_face_split(mats::AbstractMatrix...)

Generalized row-wise unique face-splitting for any number of matrices.
"""
function unique_face_split(mats::AbstractMatrix...)
    L = length(mats)
    L ≥ 1 || throw(ArgumentError("need at least one matrix"))
    nr = size(mats[1], 1)
    for M in mats
        size(M, 1) == nr ||
          throw(ArgumentError("all matrices must have the same number of rows"))
    end
    if L == 1
        return unique_face_split(mats[1])
    elseif L == 2
        return unique_face_split(mats[1], mats[2])
    else
        rows = [ reduce(unique_kronecker, (M[j, :] for M in mats)) for j in 1:nr ]
        return vcat((transpose(r) for r in rows)...)
    end
end

"""
    unique_face_split(A::AbstractMatrix, d::Integer)

Raise each row of A to the unique-Kronecker power d, stacking the
transposed results.
"""
function unique_face_split(A::AbstractMatrix, d::Integer)
    d ≥ 1 || throw(ArgumentError("d must be at least 1"))
    nr, _ = size(A)
    rows = [ transpose(unique_kronecker_power(A[j, :], d)) for j in 1:nr ]
    return vcat(rows...)
end

"""
    A ⧁ B

Operator form of unique-face-splitting product and its variants.
"""
⧁(A::AbstractMatrix, B::AbstractMatrix)      = unique_face_split(A, B)
⧁(A::AbstractMatrix)                         = unique_face_split(A)
⧁(mats::AbstractMatrix...)                   = unique_face_split(mats...)
⧁(A::AbstractMatrix, d::Integer)             = unique_face_split(A, d)