export DiffMatCache, diffmat_blocks, diffmat, unique_kronecker_jacobian, unique_kronecker_jacobian!
export unique_kronecker_jacobian2, unique_kronecker_jacobian3
export coupling_matrix


# Helper: collect combinations safely (UniqueCombinationIterator mutates in-place)
"""
    _collect_combinations(n::Int, p::Int) -> Vector{Vector{Int}}

Collect all combinations with repetition of length `p` from `{1,…,n}` in
non-decreasing order.  Each combination is an independent `Vector{Int}`.
"""
function _collect_combinations(n::Int, p::Int)
    result = Vector{Vector{Int}}(undef, binomial(n + p - 1, p))
    idx = 1
    for comb in UniqueCombinationIterator(n, p)
        result[idx] = copy(comb)
        idx += 1
    end
    return result
end


# Differentiation matrices  N_d^{(i)}  and  N^{(i)}
"""
    diffmat_blocks(n::Int, i::Int) -> Vector{SparseMatrixCSC{Int,Int}}

Compute the `n` sparse differentiation matrices ``N_d^{(i)} \\in \\mathbb{R}^{q_i \\times q_{i-1}}``
for ``d = 1, \\dots, n``, as defined in Proposition 3.6.

Each matrix satisfies ``[N_d^{(i)}]_{\\alpha,\\beta} = \\alpha_d`` when
``\\beta = \\alpha - e_d`` and ``\\alpha_d \\ge 1``, and zero otherwise.

The column `d` of the Jacobian ``\\nabla_{\\hat{s}} p_i^u(\\hat{s})`` is then
``N_d^{(i)} \\, p_{i-1}^u(\\hat{s})``.

## Arguments
- `n::Int`: dimension of the vector (number of variables).
- `i::Int`: degree of the unique monomial map (must be ≥ 1).

## Returns
- `Nblocks::Vector{SparseMatrixCSC{Int,Int}}`: vector of length `n`.

## Example
```julia-repl
julia> Nblocks = diffmat_blocks(2, 2)
2-element Vector{SparseMatrixCSC{Int64, Int64}}:
 ...
julia> Nblocks[1]   # N_1^{(2)}
3×2 SparseMatrixCSC{Int64, Int64} with 2 stored entries:
 2  ⋅
 ⋅  1
 ⋅  ⋅
julia> Nblocks[2]   # N_2^{(2)}
3×2 SparseMatrixCSC{Int64, Int64} with 2 stored entries:
 ⋅  ⋅
 1  ⋅
 ⋅  2
```
"""
function diffmat_blocks(n::Int, i::Int)
    i ≥ 1 || throw(ArgumentError("degree i must be at least 1"))
    
    qi   = binomial(n + i - 1, i)       # number of unique degree-i monomials (rows)
    qi_1 = binomial(n + i - 2, i - 1)   # number of unique degree-(i-1) monomials (cols)

    # Special case: degree 1 → each N_d^{(1)} is a q1×q0 = n×1 matrix
    # with a single 1 in row d (derivative of ŝ_d w.r.t. ŝ_d, and p_0^u = [1])
    if i == 1
        blocks = Vector{SparseMatrixCSC{Int,Int}}(undef, n)
        for d in 1:n
            blocks[d] = sparse([d], [1], [1], n, 1)
        end
        return blocks
    end

    # Build lookup: sorted combination of length (i-1) → position index
    combs_prev = _collect_combinations(n, i - 1)
    comb_to_idx = Dict{Vector{Int}, Int}()
    for (idx, c) in enumerate(combs_prev)
        comb_to_idx[c] = idx
    end

    # Preallocate COO arrays for each block
    block_rows = [Int[] for _ in 1:n]
    block_cols = [Int[] for _ in 1:n]
    block_vals = [Int[] for _ in 1:n]

    # Iterate over degree-i combinations (rows of the Jacobian)
    row = 0
    for comb in UniqueCombinationIterator(n, i)
        row += 1
        # For each variable d, check if d appears in comb
        d_prev = 0  # skip duplicate entries for same d in one comb
        for pos in eachindex(comb)
            d = comb[pos]
            d == d_prev && continue   # already processed this variable for this row
            d_prev = d

            # Count occurrences → α_d
            alpha_d = count(==(d), comb)

            # Build the reduced combination (remove one occurrence of d)
            reduced = copy(comb)
            deleteat!(reduced, pos)

            # Look up column index in I_{i-1}
            col = comb_to_idx[reduced]

            push!(block_rows[d], row)
            push!(block_cols[d], col)
            push!(block_vals[d], alpha_d)
        end
    end

    # Assemble sparse matrices
    blocks = Vector{SparseMatrixCSC{Int,Int}}(undef, n)
    for d in 1:n
        blocks[d] = sparse(block_rows[d], block_cols[d], block_vals[d], qi, qi_1)
    end
    return blocks
end


"""
    diffmat(n::Int, i::Int) -> N::SparseMatrixCSC{Int,Int}

Compute the full differentiation matrix ``N^{(i)} = [N_1^{(i)}, \\dots, N_n^{(i)}]
\\in \\mathbb{R}^{q_i \\times n \\, q_{i-1}}``, the horizontal concatenation of all
per-variable differentiation blocks.

The Jacobian of the unique monomial map satisfies (Proposition 3.6):
```math
\\nabla_{\\hat{s}} p_i^u(\\hat{s}) = N^{(i)} \\bigl(I_n \\otimes p_{i-1}^u(\\hat{s})\\bigr).
```

## Arguments
- `n::Int`: dimension of the vector.
- `i::Int`: degree of the unique monomial map (must be ≥ 1).

## Returns
- `N::SparseMatrixCSC{Int,Int}`: sparse matrix of size ``q_i \\times n q_{i-1}``.
"""
function diffmat(n::Int, i::Int)
    blocks = diffmat_blocks(n, i)
    return hcat(blocks...)
end


# DiffMatCache: precompute and store all differentiation blocks 
"""
    DiffMatCache(n::Int, pmax::Int)

Precompute and cache the sparse differentiation matrix blocks ``N_d^{(i)}`` for
``d = 1,\\dots,n`` and ``i = 1,\\dots, p_{\\max}``.  These are constant matrices that
depend only on the variable dimension `n` and the degree, so they can be assembled
once offline and reused for every Jacobian evaluation.

## Fields
- `n::Int`: variable dimension.
- `pmax::Int`: maximum polynomial degree cached.
- `blocks::Vector{Vector{SparseMatrixCSC{Int,Int}}}`:
  `blocks[i]` is a length-`n` vector of the `N_d^{(i)}` blocks for degree `i`.

## Example
```julia-repl
julia> cache = DiffMatCache(3, 4)   # n=3 variables, up to degree 4
julia> cache.blocks[2]              # the three N_d^{(2)} blocks
julia> cache.blocks[2][1]           # N_1^{(2)}
```
"""
struct DiffMatCache
    n::Int
    pmax::Int
    blocks::Vector{Vector{SparseMatrixCSC{Int,Int}}}    # blocks[i][d] = N_d^{(i)}

    function DiffMatCache(n::Int, pmax::Int)
        n ≥ 1 || throw(ArgumentError("n must be at least 1"))
        pmax ≥ 1 || throw(ArgumentError("pmax must be at least 1"))
        blks = [diffmat_blocks(n, i) for i in 1:pmax]
        return new(n, pmax, blks)
    end
end

Base.show(io::IO, c::DiffMatCache) =
    print(io, "DiffMatCache(n=$(c.n), pmax=$(c.pmax))")


# Jacobian evaluation: ∇_ŝ pᵤⁱ(ŝ)
"""
    unique_kronecker_jacobian(x::AbstractVector{T}, i::Int,
                              [cache::DiffMatCache]) -> J::Matrix{T}

Compute the Jacobian ``\\nabla_{\\hat{s}} p_i^u(\\hat{s})`` of the unique
degree-`i` monomial map evaluated at `x`, returning a dense ``q_i \\times n``
matrix.

The computation uses the formula from Proposition 3.6:
```math
\\nabla_{\\hat{s}} p_i^u(\\hat{s}) = N^{(i)} \\bigl(I_n \\otimes p_{i-1}^u(\\hat{s})\\bigr),
```
evaluated column-by-column as ``N_d^{(i)} \\, p_{i-1}^u(\\hat{s})`` for ``d = 1,\\dots,n``
to avoid forming the full Kronecker product.

## Arguments
- `x::AbstractVector{T}`: the point at which to evaluate (length `n`).
- `i::Int`: degree of the monomial map (must be ≥ 1).
- `cache::DiffMatCache` (optional): precomputed differentiation matrices.
  If omitted, they are computed on the fly.

## Returns
- `J::Matrix{T}`: the ``q_i \\times n`` Jacobian matrix.

## Example
```julia-repl
julia> x = [3.0, 5.0]
julia> J = unique_kronecker_jacobian(x, 2)
3×2 Matrix{Float64}:
 6.0  0.0
 5.0  3.0
 0.0  10.0
```
"""
function unique_kronecker_jacobian(x::AbstractVector{T}, i::Int) where {T}
    n = length(x)
    blocks = diffmat_blocks(n, i)
    return _eval_jacobian(x, i, blocks)
end

function unique_kronecker_jacobian(x::AbstractVector{T}, i::Int, cache::DiffMatCache) where {T}
    n = length(x)
    n == cache.n || throw(DimensionMismatch("vector length $n ≠ cache dimension $(cache.n)"))
    i ≤ cache.pmax || throw(ArgumentError("degree $i exceeds cache maximum $(cache.pmax)"))
    return _eval_jacobian(x, i, cache.blocks[i])
end


# Internal: column-by-column evaluation
"""
    _eval_jacobian(x, i, Nblocks) -> Matrix

Evaluate ∇_ŝ pᵢᵘ(ŝ) at x using precomputed differentiation blocks.
Column d = Nblocks[d] * p_{i-1}^u(x).  Avoids forming I_n ⊗ p_{i-1}^u(x).
"""
function _eval_jacobian(x::AbstractVector{T}, i::Int, Nblocks::Vector{<:SparseMatrixCSC}) where {T}
    n = length(x)
    qi = binomial(n + i - 1, i)

    # Degree-1 special case: Jacobian is just I_n
    if i == 1
        return Matrix{T}(I, n, n)
    end

    # Evaluate the (i-1)-th unique monomial vector
    pu_prev = _unique_monomial_vec(x, i - 1)

    # Build the Jacobian column by column: J[:, d] = N_d^{(i)} * p_{i-1}^u(x)
    J = Matrix{T}(undef, qi, n)
    @inbounds for d in 1:n
        mul!(view(J, :, d), Nblocks[d], pu_prev)
    end
    return J
end


"""
    _unique_monomial_vec(x::AbstractVector{T}, p::Int) -> Vector{T}

Evaluate the unique monomial vector ``p_p^u(x)`` of degree `p`.  Dispatches
to the specialized routines for p ≤ 4, and falls back to `unique_kronecker_power`
otherwise.  Returns a plain `Vector{T}` in all cases.
"""
@inline function _unique_monomial_vec(x::AbstractVector{T}, p::Int) where {T}
    if p == 0
        return T[one(T)]
    elseif p == 1
        return collect(x)
    elseif p == 2
        return unique_kronecker(x, x)
    elseif p == 3
        return unique_kronecker(x, x, x)
    elseif p == 4
        return unique_kronecker(x, x, x, x)
    else
        return unique_kronecker_power(x, p)
    end
end


# In-place Jacobian evaluation
"""
    unique_kronecker_jacobian!(J::AbstractMatrix{T}, x::AbstractVector{T}, i::Int,
                               cache::DiffMatCache) -> J

In-place version of [`unique_kronecker_jacobian`](@ref).  Writes into the
pre-allocated ``q_i \\times n`` matrix `J`.

## Arguments
- `J::AbstractMatrix{T}`: output matrix of size ``q_i \\times n``.
- `x::AbstractVector{T}`: evaluation point (length `n`).
- `i::Int`: degree of the monomial map.
- `cache::DiffMatCache`: precomputed differentiation matrices.

## Returns
- `J` (modified in place).
"""
function unique_kronecker_jacobian!(J::AbstractMatrix{T}, x::AbstractVector{T}, i::Int,
                                    cache::DiffMatCache) where {T}
    n = length(x)
    n == cache.n || throw(DimensionMismatch("vector length $n ≠ cache dimension $(cache.n)"))
    i ≤ cache.pmax || throw(ArgumentError("degree $i exceeds cache maximum $(cache.pmax)"))
    qi = binomial(n + i - 1, i)
    size(J) == (qi, n) || throw(DimensionMismatch("J must be $qi × $n, got $(size(J))"))

    if i == 1
        J .= zero(T)
        @inbounds for d in 1:n
            J[d, d] = one(T)
        end
        return J
    end

    pu_prev = _unique_monomial_vec(x, i - 1)
    Nblocks = cache.blocks[i]
    @inbounds for d in 1:n
        mul!(view(J, :, d), Nblocks[d], pu_prev)
    end
    return J
end


# Specialized fast Jacobians for degrees 2, 3, 4 (no sparse matmul needed)
"""
    unique_kronecker_jacobian2(x::AbstractVector{T}) -> Matrix{T}

Specialized Jacobian of ``p_2^u(\\hat{s}) = [\\hat{s}_1^2, \\hat{s}_1\\hat{s}_2, \\dots]``
for degree 2, computed directly without differentiation matrices.
"""
function unique_kronecker_jacobian2(x::AbstractVector{T}) where {T}
    n = length(x)
    q2 = n * (n + 1) ÷ 2
    J = zeros(T, q2, n)
    row = 1
    @inbounds for a in 1:n
        for b in a:n
            if a == b
                J[row, a] = 2 * x[a]
            else
                J[row, a] = x[b]
                J[row, b] = x[a]
            end
            row += 1
        end
    end
    return J
end

"""
    unique_kronecker_jacobian3(x::AbstractVector{T}) -> Matrix{T}

Specialized Jacobian of ``p_3^u(\\hat{s})`` for degree 3, computed directly.
"""
function unique_kronecker_jacobian3(x::AbstractVector{T}) where {T}
    n = length(x)
    q3 = binomial(n + 2, 3)
    J = zeros(T, q3, n)
    row = 1
    @inbounds for a in 1:n
        for b in a:n
            for c in b:n
                # Monomial x[a]*x[b]*x[c], derivatives via product rule
                # Count multiplicities for each unique index
                if a == b == c
                    J[row, a] = 3 * x[a] * x[a]
                elseif a == b  # a == b ≠ c
                    J[row, a] = 2 * x[a] * x[c]
                    J[row, c] = x[a] * x[a]
                elseif b == c  # a ≠ b == c
                    J[row, a] = x[b] * x[b]
                    J[row, b] = 2 * x[a] * x[b]
                elseif a == c  # shouldn't happen for a ≤ b ≤ c with a ≠ b, but just in case
                    J[row, a] = 2 * x[a] * x[b]
                    J[row, b] = x[a] * x[a]
                else  # all distinct
                    J[row, a] = x[b] * x[c]
                    J[row, b] = x[a] * x[c]
                    J[row, c] = x[a] * x[b]
                end
                row += 1
            end
        end
    end
    return J
end


# Matrix (snapshot) interface: Jacobian at each column
"""
    unique_kronecker_jacobian(X::AbstractMatrix{T}, i::Int,
                              [cache::DiffMatCache]) -> Vector{Matrix{T}}

Compute the Jacobian of the unique degree-`i` monomial map at each column of `X`,
returning a vector of ``q_i \\times n`` Jacobian matrices.

## Arguments
- `X::AbstractMatrix{T}`: matrix whose columns are the evaluation points (size ``n \\times k``).
- `i::Int`: degree of the monomial map.
- `cache::DiffMatCache` (optional): precomputed differentiation matrices.

## Returns
- `Js::Vector{Matrix{T}}`: `Js[ℓ]` is the Jacobian at `X[:, ℓ]`.
"""
function unique_kronecker_jacobian(X::AbstractMatrix{T}, i::Int) where {T}
    n = size(X, 1)
    cache = DiffMatCache(n, i)
    return unique_kronecker_jacobian(X, i, cache)
end

function unique_kronecker_jacobian(X::AbstractMatrix{T}, i::Int, cache::DiffMatCache) where {T}
    n, k = size(X)
    n == cache.n || throw(DimensionMismatch("matrix row count $n ≠ cache dimension $(cache.n)"))
    i ≤ cache.pmax || throw(ArgumentError("degree $i exceeds cache maximum $(cache.pmax)"))
    qi = binomial(n + i - 1, i)
    Js = [Matrix{T}(undef, qi, n) for _ in 1:k]
    @inbounds for ℓ in 1:k
        unique_kronecker_jacobian!(Js[ℓ], view(X, :, ℓ), i, cache)
    end
    return Js
end


# Coupling matrix  C(ŝ) = Σ_{i=2}^{p} G^{(i)} ∇_ŝ pᵢᵘ(ŝ)
"""
    coupling_matrix(x::AbstractVector{T}, G::Vector{<:AbstractMatrix},
                    cache::DiffMatCache; pstart::Int=2) -> C::Matrix{T}

Compute the coupling matrix
```math
C(\\hat{s}) = \\sum_{i=p_{\\text{start}}}^{p} G^{(i)} \\nabla_{\\hat{s}} p_i^u(\\hat{s})
\\in \\mathbb{R}^{j \\times n}
```
where ``G^{(i)} = V_j^\\top V^{(i)} \\in \\mathbb{R}^{j \\times q_i}`` are precomputed
projection matrices and `pstart` (default 2) is the starting degree.

## Arguments
- `x::AbstractVector{T}`: reduced state ``\\hat{s}`` (length `n`).
- `G::Vector{<:AbstractMatrix}`: `G[i-pstart+1]` is ``G^{(i)}`` of size ``j \\times q_i``.
- `cache::DiffMatCache`: precomputed differentiation matrices.
- `pstart::Int`: first degree in the summation (default 2).

## Returns
- `C::Matrix{T}`: the ``j \\times n`` coupling matrix.
"""
function coupling_matrix(x::AbstractVector{T}, G::Vector{<:AbstractMatrix},
                         cache::DiffMatCache; pstart::Int=2) where {T}
    n = length(x)
    n == cache.n || throw(DimensionMismatch("vector length $n ≠ cache dimension $(cache.n)"))
    length(G) ≥ 1 || throw(ArgumentError("need at least one projection matrix G"))

    j = size(G[1], 1)
    p = pstart + length(G) - 1
    p ≤ cache.pmax || throw(ArgumentError("degree $p exceeds cache maximum $(cache.pmax)"))

    C = zeros(T, j, n)
    Ji = Matrix{T}(undef, 0, 0)  # will be resized
    for i in pstart:p
        qi = binomial(n + i - 1, i)
        if size(Ji) != (qi, n)
            Ji = Matrix{T}(undef, qi, n)
        end
        unique_kronecker_jacobian!(Ji, x, i, cache)
        # C += G^{(i)} * J_i
        Gi = G[i - pstart + 1]
        mul!(C, Gi, Ji, one(T), one(T))   # C = 1*C + 1*Gi*Ji
    end
    return C
end

"""
    coupling_matrix(X::AbstractMatrix{T}, G::Vector{<:AbstractMatrix},
                    cache::DiffMatCache; pstart::Int=2) -> Vector{Matrix{T}}

Batch version: compute the coupling matrix at each column of `X`.
"""
function coupling_matrix(X::AbstractMatrix{T}, G::Vector{<:AbstractMatrix},
                         cache::DiffMatCache; pstart::Int=2) where {T}
    k = size(X, 2)
    return [coupling_matrix(view(X, :, ℓ), G, cache; pstart=pstart) for ℓ in 1:k]
end