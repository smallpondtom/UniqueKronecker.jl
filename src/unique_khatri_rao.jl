export unique_khatri_rao, ⦼

"""
    unique_khatri_rao(A::AbstractMatrix, B::AbstractMatrix)

Column-wise unique Kronecker (Khatri-Rao) of A and B.
"""
function unique_khatri_rao(A::AbstractMatrix, B::AbstractMatrix)
    size(A,2) == size(B,2) ||
      throw(ArgumentError("matrices must have same number of columns"))
    return hcat(map(unique_kronecker, eachcol(A), eachcol(B))...)
end

unique_khatri_rao(A::AbstractMatrix) = unique_khatri_rao(A, A)

"""
    unique_khatri_rao(mats::AbstractMatrix...)

Generalized column-wise unique Kronecker of any number of matrices.
"""
function unique_khatri_rao(mats::AbstractMatrix...)
    L = length(mats)
    L ≥ 1 || throw(ArgumentError("need at least one matrix"))
    nc = size(mats[1], 2)
    for M in mats
        size(M,2) == nc ||
          throw(ArgumentError("all matrices must have the same number of columns"))
    end
    if L == 1
        return unique_khatri_rao(mats[1])
    elseif L == 2
        return unique_khatri_rao(mats[1], mats[2])
    else
        cols = [ reduce(unique_kronecker, (M[:,j] for M in mats)) for j in 1:nc ]
        return hcat(cols...)
    end
end

"""
    unique_khatri_rao(A::AbstractMatrix, d::Integer)

Raise each column of A to the unique-Kronecker power d.
"""
function unique_khatri_rao(A::AbstractMatrix, d::Integer)
    d ≥ 1 || throw(ArgumentError("d must be at least 1"))
    nc = size(A,2)
    cols = [ unique_kronecker_power(A[:,j], d) for j in 1:nc ]
    return hcat(cols...)
end

# operator aliases
⦼(A::AbstractMatrix, B::AbstractMatrix)    = unique_khatri_rao(A, B)
⦼(mats::AbstractMatrix...)                 = unique_khatri_rao(mats...)
⦼(A::AbstractMatrix, d::Integer)           = unique_khatri_rao(A, d)