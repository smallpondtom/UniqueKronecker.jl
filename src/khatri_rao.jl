export khatri_rao, ⊙

"""
    khatri_rao(A::AbstractMatrix, B::AbstractMatrix)

Compute the Khatri-Rao (column-wise Kronecker) product of matrices `A` and `B`.

Given matrices `A` (m×n) and `B` (p×n) with the same number of columns,
returns an (mp×n) matrix where each column is the Kronecker product of the
corresponding columns of `A` and `B`.

# Arguments
- `A::AbstractMatrix`: First input matrix
- `B::AbstractMatrix`: Second input matrix (must have same number of columns as `A`)

# Returns
- Matrix of size (m*p × n) containing column-wise Kronecker products

# Examples
```julia
A = [1 2; 3 4]
B = [5 6; 7 8]
kr = khatri_rao(A, B)
# Equivalent to: hcat(kron(A[:,1], B[:,1]), kron(A[:,2], B[:,2]))
```

See also: [`unique_khatri_rao`](@ref), [`⊙`](@ref)
"""
function khatri_rao(A::AbstractMatrix, B::AbstractMatrix)
    size(A,2) == size(B,2) ||
      throw(ArgumentError("matrices must have same number of columns"))
    return hcat(map(kronecker, eachcol(A), eachcol(B))...)
end

"""
    khatri_rao(A::AbstractMatrix)

Compute the Khatri-Rao product of matrix `A` with itself.

# Examples
```julia
A = [1 2; 3 4]
kr = khatri_rao(A)  # Same as khatri_rao(A, A)
```
"""
khatri_rao(A::AbstractMatrix) = khatri_rao(A, A)

"""
    ⊙(A::AbstractMatrix, B::AbstractMatrix)
    ⊙(A::AbstractMatrix)

Unicode operator for the Khatri-Rao product. Type `\\odot<TAB>` to enter.

# Examples
```julia
A = [1 2; 3 4]
B = [5 6; 7 8]
kr = A ⊙ B  # Same as khatri_rao(A, B)
```

See also: [`khatri_rao`](@ref)
"""
⊙(A::AbstractMatrix, B::AbstractMatrix) = khatri_rao(A, B)
⊙(A::AbstractMatrix) = khatri_rao(A, A)

"""
    khatri_rao(mats::AbstractMatrix...)

Generalized Khatri-Rao product for multiple matrices.

Computes the column-wise Kronecker product across all input matrices.
All matrices must have the same number of columns.

# Arguments
- `mats::AbstractMatrix...`: Variable number of matrices with matching column counts

# Examples
```julia
A = [1 2; 3 4]
B = [5 6; 7 8]
C = [9 10; 11 12]
kr = khatri_rao(A, B, C)
kr = A ⊙ B ⊙ C  # Using operator
```
"""
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

"""
    khatri_rao(A::AbstractMatrix, d::Integer)

Compute the `d`-th Khatri-Rao power of matrix `A`.

Each column of the result is the `d`-th Kronecker power of the corresponding
column of `A`.

# Arguments
- `A::AbstractMatrix`: Input matrix
- `d::Integer`: Power (must be at least 1)

# Examples
```julia
A = [1 2; 3 4]
kr2 = khatri_rao(A, 2)  # Each column: kron(A[:,j], A[:,j])
kr2 = A ⊙ 2             # Using operator
```
"""
function khatri_rao(A::AbstractMatrix, d::Integer)
    d ≥ 1 || throw(ArgumentError("d must be at least 1"))
    n = size(A, 2)
    cols = [ reduce(kron, ntuple(_->A[:, j], d)) for j in 1:n ]
    return hcat(cols...)
end

⊙(A::AbstractMatrix, d::Integer) = khatri_rao(A, d)