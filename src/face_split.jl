export face_split, ⊖

"""
    face_split(A::AbstractMatrix, B::AbstractMatrix)

Row-wise "face-splitting" product of A and B.  Each pair of rows
A[j,:], B[j,:] is Kronenecker-multiplied and transposed, then
all results are stacked vertically.
"""
function face_split(A::AbstractMatrix, B::AbstractMatrix)
    size(A,1) == size(B,1) ||
      throw(ArgumentError("matrices must have same number of rows"))
    return vcat(map(transpose ∘ kron, eachrow(A), eachrow(B))...)
end

face_split(A::AbstractMatrix) = face_split(A, A)

"""
    face_split(mats::AbstractMatrix...)

Generalized row-wise face_split of any number of matrices.
"""
function face_split(mats::AbstractMatrix...)
    L = length(mats)
    L ≥ 1 || throw(ArgumentError("need at least one matrix"))
    nr = size(mats[1], 1)
    for M in mats
        size(M,1) == nr ||
          throw(ArgumentError("all matrices must have the same number of rows"))
    end
    if L == 1
        return face_split(mats[1])
    elseif L == 2
        return face_split(mats[1], mats[2])
    else
        rows = [ reduce(kron, (M[j, :] for M in mats)) for j in 1:nr ]
        return vcat((transpose(r) for r in rows)...)
    end
end

"""
    face_split(A::AbstractMatrix, d::Integer)

Raise each row of A to the face-split (Kronecker) power d.
"""
function face_split(A::AbstractMatrix, d::Integer)
    d ≥ 1 || throw(ArgumentError("d must be at least 1"))
    nr, _ = size(A)
    rows = [ transpose(reduce(kron, ntuple(_->A[j, :], d))) for j in 1:nr ]
    return vcat(rows...)
end

# operator aliases
⊖(A::AbstractMatrix, B::AbstractMatrix)      = face_split(A, B)
⊖(A::AbstractMatrix)                         = face_split(A)
⊖(mats::AbstractMatrix...)                   = face_split(mats...)
⊖(A::AbstractMatrix, d::Integer)             = face_split(A, d)

