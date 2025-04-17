export khatri_rao, ⊙

function khatri_rao(A::AbstractMatrix, B::AbstractMatrix)
    size(A,2) == size(B,2) ||
      throw(ArgumentError("matrices must have same number of columns"))
    return hcat(map(kronecker, eachcol(A), eachcol(B))...)
end

khatri_rao(A::AbstractMatrix) = khatri_rao(A, A)

⊙(A::AbstractMatrix, B::AbstractMatrix) = khatri_rao(A, B)
⊙(A::AbstractMatrix) = khatri_rao(A, A)

function khatri_rao(mats::AbstractMatrix...)
    L = length(mats)
    L ≥ 1 || throw(ArgumentError("need at least one matrix"))
    nc = size(mats[1], 2)
    for M in mats
        size(M, 2) == nc || 
            throw(ArgumentError("all matrices must have the same number of columns"))
    end
    if L == 1
        return khatri_rao(mats[1])
    elseif L == 2
        return khatri_rao(mats[1], mats[2])
    else
        cols = [ reduce(kronecker, (M[:, j] for M in mats)) for j in 1:nc ]
    end
    return hcat(cols...)
end

⊙(mats::AbstractMatrix...) = khatri_rao(mats...)

function khatri_rao(A::AbstractMatrix, d::Integer)
    d ≥ 1 || throw(ArgumentError("d must be at least 1"))
    n = size(A, 2)
    cols = [ reduce(kron, ntuple(_->A[:, j], d)) for j in 1:n ]
    return hcat(cols...)
end

⊙(A::AbstractMatrix, d::Integer) = khatri_rao(A, d)