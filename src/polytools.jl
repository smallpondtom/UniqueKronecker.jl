export dupmat, symmtzrmat, elimat, commat, eliminate, duplicate, duplicate_symmetric
export make_poly_op
export kron_snapshot_matrix, unique_kron_snapshot_matrix

"""
    dupmat(n::Int, p::Int) -> Dp::SparseMatrixCSC{Int}

Create a duplication matrix of order `p` for a vector of length `n` [MagnusEML1980](@citet).

## Arguments
- `n::Int`: The length of the vector.
- `p::Int`: The order of the duplication matrix, e.g., `p = 2` for x ⊗ x.

## Output
- `Dp::SparseMatrixCSC{Int}`: The duplication matrix of order `p`.

## Example
```julia-repl
julia> dupmat(2,2)
4×3 SparseMatrixCSC{Int64, Int64} with 4 stored entries:
 1  ⋅  ⋅
 ⋅  1  ⋅
 ⋅  1  ⋅
 ⋅  ⋅  1
```
"""
function dupmat(n::Int, p::Int)
    if p == 2  # fast algorithm for quadratic case
        m = n * (n + 1) / 2
        nsq = n^2
        r = 1
        a = 1
        v = zeros(nsq)
        cn = cumsum(n:-1:2)
        for i = 1:n
            v[r:(r+i-2)] = (i - n) .+ cn[1:(i-1)]
            r = r + i - 1

            v[r:r+n-i] = a:a.+(n-i)
            r = r + n - i + 1
            a = a + n - i + 1
        end
        D = sparse(1:nsq, v, ones(length(v)), nsq, m)
        return D
    else  # general algorithm for any order p
        # Calculate the number of unique elements in a symmetric tensor of order p
        num_unique_elements = binomial(n + p - 1, p)
        
        # Create the duplication matrix with appropriate size
        Dp = spzeros(Int, n^p, num_unique_elements)

        function elements!(D, l, indices)
            perms = unique([sum((indices[σ[i]] - 1) * n^(p - i) for i in 1:p) + 1 for σ in permutations(1:p)])
            
            @inbounds for p in perms
                D[p, l] = 1
            end
        end
        elements!(D, l, indices...) = elements!(D, l, indices)

        # function generate_nonredundant_combinations(n::Int, p::Int)
        #     if p == 1
        #         return [[i] for i in 1:n]
        #     else
        #         lower_combs = generate_nonredundant_combinations(n, p - 1)
        #         return [vcat(comb, [i]) for comb in lower_combs for i in comb[end]:n]
        #     end
        # end

        # Generate all combinations (i1, i2, ..., ip) with i1 ≤ i2 ≤ ... ≤ ip
        # combs = generate_nonredundant_combinations(n, p)
        combs = with_replacement_combinations(1:n, p) 

        # INFO: Vectorization involves overhead of `reduce(hcat, combs)`
        # combs = reduce(hcat, combs)
        # combs = Vector{eltype(combs)}[eachrow(combs)...]
        # # Fill in the duplication matrix using the computed combinations
        # elements!.(Ref(Dp), 1:num_unique_elements, combs...)

        # Remove overhead of `reduce(hcat, combs)` by using a loop
        l = 1
        for indices in combs
            elements!(Dp, l, indices)
            l += 1
        end

        return Dp
    end
end


"""
    symmtzrmat(n::Int, p::Int) -> Sp::SparseMatrixCSC{Float64}

Create a symmetrizer matrix of order `p` for a vector of length `n` [MagnusEML1980](@citet).

## Arguments
- `n::Int`: The length of the vector.
- `p::Int`: The order of the symmetrizer matrix, e.g., `p = 2` for x ⊗ x.

## Output
- `Sp::SparseMatrixCSC{Float64}`: The symmetrizer matrix of order `p`.

## Example
```julia-repl
julia> symmtzrmat(2,2)
4×4 SparseMatrixCSC{Float64, Int64} with 6 stored entries:
 1.0   ⋅    ⋅    ⋅ 
  ⋅   0.5  0.5   ⋅
  ⋅   0.5  0.5   ⋅
  ⋅    ⋅    ⋅   1.0
```
"""
function symmtzrmat(n::Int, p::Int)
    if p == 2  # fast algorithm for quadratic case
        np = n^p
        return 0.5 * (sparse(1.0I, np, np) + commat(n, n))
    else
        # Create the symmetrizer matrix with appropriate size
        np = Int(n^p)
        Sp = spzeros(Float64, np, np)

        function elements!(N, l, indices)
            perms = [sum((indices[σ[i]] - 1) * n^(p - i) for i in 1:p) + 1 for σ in permutations(1:p)]
            
            # For cases where two or all indices are the same, 
            # we should not count permutations more than once.
            unique_perms = countmap(perms)

            # Assign the column to the matrix N
            @inbounds for (perm, count) in unique_perms
                N[perm, l] = count / factorial(p)
            end
        end
        elements!(N, l, indices...) = elements!(N, l, indices)

        function generate_redundant_combinations(n::Int, p::Int)
            iterators = ntuple(_ -> 1:n, p)
            return [collect(product) for product in Iterators.product(iterators...)]
        end

        # Generate all combinations (i1, i2, ..., ip) with i1 ≤ i2 ≤ ... ≤ ip
        combs = generate_redundant_combinations(n, p)

        # INFO: Vectorization involves overhead of `reduce(hcat, combs)`
        # combs = reduce(hcat, combs)
        # combs = Vector{eltype(combs)}[eachrow(combs)...]
        # # Fill in the duplication matrix using the computed combinations
        # elements!.(Ref(Sp), 1:np, combs...)

        # Remove overhead of `reduce(hcat, combs)` by using a loop
        l = 1
        for indices in combs
            elements!(Sp, l, indices)
            l += 1
        end
        
        return Sp
    end
end

# INFO: EXPERIMENTAL (THEORETICALLY INCORRECT)
# function symmtzrmat2(n::Int, p::Int)
#     np = n^p
#     Sp = sparse(1.0I, np, np)

#     for i in 1:p-1
#         j = p - i
#         Sp .+= commat(n^i, n^j)
#     end
#     return Sp / p
# end


"""
    elimat(n::Int, p::Int) -> Lp::SparseMatrixCSC{Int}

Create an elimination matrix of order `p` for a vector of length `n` [MagnusEML1980](@citet).

## Arguments
- `n::Int`: The length of the vector.
- `p::Int`: The order of the elimination matrix, e.g., `p = 2` for x ⊗ x.

## Output
- `Lp::SparseMatrixCSC{Int}`: The elimination matrix of order `p`.

## Example
```julia-repl
julia> elimat(2,2)
3×4 SparseMatrixCSC{Int64, Int64} with 3 stored entries:
 1  ⋅  ⋅  ⋅
 ⋅  1  ⋅  ⋅
 ⋅  ⋅  ⋅  1
```
"""
function elimat(n::Int, p::Int)
    if p == 2  # fast algorithm for quadratic case
        T = tril(ones(n, n)) # Lower triangle of 1's
        f = findall(x -> x == 1, T[:]) # Get linear indexes of 1's
        k = n * (n + 1) / 2 # Row size of L
        n2 = n * n # Colunm size of L
        x = f + n2 * (0:k-1) # Linear indexes of the 1's within L'

        row = [mod(a, n2) != 0 ? mod(a, n2) : n2 for a in x]
        col = [mod(a, n2) != 0 ? div(a, n2) + 1 : div(a, n2) for a in x]
        L = sparse(row, col, ones(length(x)), n2, k)
        L = L' # Now transpose to actual L
        return L
    else   # general algorithm for any order p
        # Calculate the number of rows in L
        num_rows = binomial(n + p - 1, p)

        # Initialize the output matrix L
        Lp = spzeros(Int, num_rows, n^p)

        # Generate all combinations with repetition of n elements from {1, 2, ..., n}
        combs = with_replacement_combinations(1:n, p) 

        # Fill the matrix L
        for (l, comb) in enumerate(combs)
            v = [1]  # Start with a scalar 1
            @inbounds for d in 1:p
                e = collect(1:n .== comb[d])  # Create the indicator vector
                v = v ⊗ e  # Build the Kronecker product
            end
            Lp[l, :] = v  # Assign the row
        end

        return Lp
    end
end


"""
    commat(m::Int, n::Int) → K

Create commutation matrix `K` of dimension `m x n` [MagnusEML1980](@citet).

## Arguments
- `m::Int`: row dimension of the commutation matrix
- `n::Int`: column dimension of the commutation matrix

## Returns
- `K`: commutation matrix

## Example
```julia-repl
julia> commat(2,2)
4×4 SparseMatrixCSC{Float64, Int64} with 4 stored entries:
 1.0   ⋅    ⋅    ⋅
  ⋅    ⋅   1.0   ⋅
  ⋅   1.0   ⋅    ⋅ 
  ⋅    ⋅    ⋅   1.0
```
"""
function commat(m::Int, n::Int)
    mn = Int(m * n)
    A = reshape(1:mn, m, n)
    v = vec(A')
    K = sparse(1.0I, mn, mn)
    K = K[v, :]
    return K
end


"""
    commat(m::Int) → K

Dispatch for the commutation matrix of dimensions (m, m)

## Arguments
- `m::Int`: row and column dimension of the commutation matrix

## Returns
- `K`: commutation matrix
"""
commat(m::Int) = commat(m, m)  # dispatch


"""
    eliminate(A::AbstractArray, p::Int)

Eliminate the redundant polynomial coefficients in the matrix `A` and return the matrix 
with unique coefficients.

## Arguments
- `A::AbstractArray`: A matrix
- `p::Int`: The order of the polynomial, e.g., `p = 2` for x ⊗ x.

## Returns
- matrix with unique coefficients

## Example
```julia-repl
julia> n = 2; P = rand(n,n); P *= P'; p = vec(P)
4-element Vector{Float64}:
 0.5085988756090203
 0.7704767970682769
 0.7704767970682769
 1.310279680309927

julia> Q = rand(n,n); Q *= Q'; q = vec(Q)
4-element Vector{Float64}:
 0.40940214810208353
 0.2295272821417254
 0.2295272821417254
 0.25503767587483905

julia> A2 = [p'; q']
2×4 Matrix{Float64}:
 0.257622  0.44721   0.0202203  0.247649
 1.55077   0.871029  0.958499   0.650717

julia> eliminate(A2, 2)
2×3 Matrix{Float64}:
 0.508599  1.54095   1.31028
 0.409402  0.459055  0.255038
```
"""
function eliminate(A::AbstractArray, p::Int)
    n = size(A, 1)
    Dp = dupmat(n, p)
    return A * Dp
end


"""
    duplicate(A::AbstractArray)

Duplicate the redundant polynomial coefficients in the matrix `A` with a unique set of coefficients
and return the matrix with redundant coefficients.

## Arguments
- `A::AbstractArray`: A matrix
- `p::Int`: The order of the polynomial, e.g., `p = 2` for x ⊗ x.

## Returns
- matrix with redundant coefficients

## Example
```julia-repl
julia> n = 2; P = rand(n,n); P *= P'; p = vec(P)
4-element Vector{Float64}:
 0.5085988756090203
 0.7704767970682769
 0.7704767970682769
 1.310279680309927

julia> Q = rand(n,n); Q *= Q'; q = vec(Q)
4-element Vector{Float64}:
 0.40940214810208353
 0.2295272821417254
 0.2295272821417254
 0.25503767587483905

julia> A2 = [p'; q']
2×4 Matrix{Float64}:
 0.257622  0.44721   0.0202203  0.247649
 1.55077   0.871029  0.958499   0.650717

julia> D2 = dupmat(2,2)
4×3 SparseMatrixCSC{Int64, Int64} with 4 stored entries:
 1  ⋅  ⋅
 ⋅  1  ⋅
 ⋅  1  ⋅
 ⋅  ⋅  1

julia> A2 * D2
2×3 Matrix{Float64}:
 0.508599  1.54095   1.31028
 0.409402  0.459055  0.255038

julia> duplicate(A2 * D2, 2)
2×4 Matrix{Float64}:
 0.508599  1.54095   0.0  1.31028
 0.409402  0.459055  0.0  0.255038
```

"""
function duplicate(A::AbstractArray, p::Int)
    n = size(A, 1)
    Lp = elimat(n, p)
    return A * Lp
end


"""
    duplicate_symmetric(A::AbstractArray, p::Int)

Duplicate the redundant polynomial coefficients in the matrix `A` with a unique set of coefficients
and return the matrix with redundant coefficients which are duplicated symmetrically.
This guarantees that the operator is symmetric. The difference from `duplicate` is that
we use the elimination matrix `Lp` and the symmetric commutation matrix `Sp` to multiply the `A` matrix.

## Arguments
- `A::AbstractArray`: A matrix
- `p::Int`: The order of the polynomial, e.g., `p = 2` for x ⊗ x.

## Returns
- matrix with redundant coefficients duplicated symmetrically

## Example
```julia-repl
julia> n = 2; P = rand(n,n); P *= P'; p = vec(P)
4-element Vector{Float64}:
 0.5085988756090203
 0.7704767970682769
 0.7704767970682769
 1.310279680309927

julia> Q = rand(n,n); Q *= Q'; q = vec(Q)
4-element Vector{Float64}:
 0.40940214810208353
 0.2295272821417254
 0.2295272821417254
 0.25503767587483905

julia> A2 = [p'; q']
2×4 Matrix{Float64}:
 0.257622  0.44721   0.0202203  0.247649
 1.55077   0.871029  0.958499   0.650717

julia> D2 = dupmat(2,2)
4×3 SparseMatrixCSC{Int64, Int64} with 4 stored entries:
 1  ⋅  ⋅
 ⋅  1  ⋅
 ⋅  1  ⋅
 ⋅  ⋅  1

julia> A2 * D2
2×3 Matrix{Float64}:
 0.508599  1.54095   1.31028
 0.409402  0.459055  0.255038

julia> duplicate_symmetric(A2 * D2, 2)
2×4 Matrix{Float64}:
 0.508599  0.770477  0.770477  1.31028
 0.409402  0.229527  0.229527  0.255038
```
"""
function duplicate_symmetric(A::AbstractArray, p::Int)
    n = size(A, 1)
    Lp = elimat(n, p)
    Sp = symmtzrmat(n, p)
    return A * Lp * Sp
end


"""
    kron_snapshot_matrix(Xmat::AbstractArray{T}, p::Int) where {T<:Number}

Take the `p`-order Kronecker product of each state of the snapshot matrix `Xmat`.

## Arguments
- `Xmat::AbstractArray{T}`: state snapshot matrix
- `p::Int`: order of the Kronecker product

## Returns
- kronecker product state snapshot matrix
"""
function kron_snapshot_matrix(Xmat::AbstractArray{T}, p::Int) where {T<:Number}
    function kron_timestep(x)
        return x[:,:] ⊗ p  # x has to be AbstractMatrix (Kronecker.jl)
    end
    tmp = kron_timestep.(eachcol(Xmat))
    return reduce(hcat, tmp)
end


"""
    unique_kron_snapshot_matrix(Xmat::AbstractArray{T}, p::Int) where {T<:Number}

Take the `p`-order unique Kronecker product of each state of the snapshot matrix `Xmat`.

## Arguments
- `Xmat::AbstractArray{T}`: state snapshot matrix
- `p::Int`: order of the Kronecker product

## Returns
- unique kronecker product state snapshot matrix
"""
function unique_kron_snapshot_matrix(Xmat::AbstractArray{T}, p::Int) where {T<:Number}
    function unique_kron_timestep(x)
        return ⊘(x, p)
    end
    tmp = unique_kron_timestep.(eachcol(Xmat))
    return reduce(hcat, tmp)
end


"""
    make_poly_op(n::Int, inds::AbstractArray{<:NTuple{P,<:Int}}, vals::AbstractArray{<:Real}; 
                    nonredundant::Bool=true, symmetric::Bool=true) where P

Helper function to construct the polynomial operator from the indices and values. The indices must
be a 1-dimensional array of, e.g., tuples of the form `(i,j)` where `i,j` are the indices of the
polynomial term. For example, for the polynomial term ``2.5x_1x_2`` for ``\\dot{x}_3`` would have an
index of `(1,2,3)` with a value of `2.5`. The `nonredundant` argument specifies which polynomial
operator to construct (the redundant or non-redundant operator). Note that the values must be a
1-dimensional array of the same length as the indices. The `symmetric` argument specifies whether
to construct the operator with symmetric coefficients.

## Arguments
- `n::Int`: dimension of the polynomial operator
- `inds::AbstractArray{<:NTuple{P,<:Int}}`: indices of the polynomial term
- `vals::AbstractArray{<:Real}`: values of the polynomial term
- `nonredundant::Bool=true`: whether to construct the non-redundant operator
- `symmetric::Bool=true`: whether to construct the symmetric operator

## Returns
- the polynomial operator
"""
function make_poly_op(n::Int, inds::AbstractArray{<:NTuple{P,<:Int}}, vals::AbstractArray{<:Real}; 
                    nonredundant::Bool=true, symmetric::Bool=true) where P
    p = P - 1
    @assert length(inds) == length(vals) "The length of indices and values must be the same."
    # Initialize a p-order tensor of size n^p
    S = zeros(ntuple(_ -> n, p+1))
    
    # nested function to assign indices and values considering symmetry
    function assign!(Smat, ind, val)
        if symmetric
            element_idx = ind[1:end-1]
            last_idx = ind[end]
            perms = unique(permutations(element_idx))
            contribution = val / length(perms)
            for perm in perms
                Smat[perm...,last_idx] = contribution
            end
        else
            Smat[ind...] = val
        end
    end

    # Vectorize the assignment process
    assign!.(Ref(S), inds, vals)

    # Flatten the p-order tensor into a matrix form with non-unique coefficients
    A = spzeros(n, n^p)

    # nested function to flatten the tensor
    function flatten!(Amat, Smat, index)
        Amat[index, :] = vec(Smat[ntuple(_ -> :, p)..., index])
    end

    # Vectorize the flattening process
    flatten!.(Ref(A), Ref(S), 1:n)

    if nonredundant
        return eliminate(A, p)
    else
        return A
    end
end

# function make_poly_op(n::Int, inds::AbstractArray{<:NTuple{P,<:Int}}, vals::AbstractArray{<:Real}; 
#                     nonredundant::Bool=true, symmetric::Bool=true) where P
#     p = P - 1
#     @assert length(inds) == length(vals) "The length of indices and values must be the same."
#     # Initialize a p-order tensor of size n^p
#     S = zeros(ntuple(_ -> n, p+1))
    
#     # Iterate over indices and assign values considering symmetry
#     for (ind, val) in zip(inds, vals)
#         if symmetric
#             element_idx = ind[1:end-1]
#             last_idx = ind[end]
#             perms = unique(permutations(element_idx))
#             contribution = val / length(perms)
#             for perm in perms
#                 S[perm...,last_idx] = contribution
#             end
#         else
#             S[ind...] = val
#         end
#     end

#     # Flatten the p-order tensor into a matrix form with non-unique coefficients
#     A = spzeros(n, n^p)
#     for i in 1:n
#         A[i, :] = vec(S[ntuple(_ -> :, p)..., i])
#     end

#     if nonredundant
#         return eliminate(A, p)
#     else
#         return A
#     end
# end


# function make_poly_op_faster(n::Int, inds::AbstractVector{<:NTuple{P, <:Int}}, vals::AbstractVector{<:Real};
#                     nonredundant::Bool=true, symmetric::Bool=true) where P
#     p = P - 1
#     @assert length(inds) == length(vals) "The length of indices and values must be the same."

#     # Initialize the tensor S as a dictionary to save memory
#     S = Dict{NTuple{P, Int}, Float64}()

#     # Precompute the LinearIndices for column indexing
#     dims = ntuple(_ -> n, p)
#     LI = LinearIndices(dims)

#     # Iterate over indices and values (sequentially to preserve order)
#     for idx in eachindex(vals)
#         ind = inds[idx]
#         val = vals[idx]
#         if symmetric
#             element_idx = ind[1:end-1]
#             last_idx = ind[end]
#             # Generate all unique permutations
#             perms = unique(permutations(element_idx))
#             contribution = val / length(perms)
#             for perm in perms
#                 key = (perm..., last_idx)
#                 # Overwrite the value to match the original function's behavior
#                 S[key] = contribution
#             end
#         else
#             key = ind
#             S[key] = val
#         end
#     end

#     # Build the sparse matrix A directly from the dictionary
#     ncols = n^p
#     A = spzeros(n, ncols)

#     for (key, val) in S
#         idxs = key
#         row = idxs[end]
#         col_idxs = idxs[1:end-1]
#         col = LI[col_idxs...]
#         A[row, col] = val  # Assign the value directly
#     end

#     if nonredundant
#         return eliminate(A, p)  # Replace with your actual implementation
#     else
#         return A
#     end
# end


# function make_poly_op_parallel(n::Int, inds::AbstractVector{<:NTuple{P, <:Int}}, vals::AbstractVector{<:Real};
#                              nonredundant::Bool=true, symmetric::Bool=true) where P
#     p = P - 1
#     @assert length(inds) == length(vals) "The length of indices and values must be the same."

#     # Number of threads
#     nthreads = Threads.nthreads()
#     if nthreads == 1
#         @warn "Only one thread is available. Using the sequential version."
#         return make_poly_op_faster(n, inds, vals, nonredundant=nonredundant, symmetric=symmetric)
#     end
#     # Initialize per-thread dictionaries to store partial results
#     S_per_thread = [Dict{NTuple{P, Int}, Tuple{Int, Float64}}() for _ in 1:nthreads]

#     # Precompute the LinearIndices for column indexing
#     dims = ntuple(_ -> n, p)
#     LI = LinearIndices(dims)

#     # Threaded loop over indices and values
#     Base.Threads.@threads for idx in eachindex(vals)
#         ind = inds[idx]
#         val = vals[idx]
#         tid = Base.Threads.threadid()
#         S_local = S_per_thread[tid]
#         if symmetric
#             element_idx = ind[1:end-1]
#             last_idx = ind[end]
#             # Generate all unique permutations
#             perms = unique(permutations(element_idx))
#             contribution = val / length(perms)
#             for perm in perms
#                 key = (perm..., last_idx)
#                 # Overwrite if idx is greater or equal
#                 if !haskey(S_local, key) || idx >= S_local[key][1]
#                     S_local[key] = (idx, contribution)
#                 end
#             end
#         else
#             key = ind
#             # Overwrite if idx is greater or equal
#             if !haskey(S_local, key) || idx >= S_local[key][1]
#                 S_local[key] = (idx, val)
#             end
#         end
#     end

#     # Merge per-thread dictionaries into a single dictionary
#     S = Dict{NTuple{P, Int}, Tuple{Int, Float64}}()
#     for S_local in S_per_thread
#         for (key, (idx, contribution)) in S_local
#             if !haskey(S, key) || idx >= S[key][1]
#                 S[key] = (idx, contribution)
#             end
#         end
#     end

#     # Build the sparse matrix A directly from the combined dictionary
#     ncols = n^p
#     A = spzeros(n, ncols)

#     for (key, (_idx, val)) in S
#         idxs = key
#         row = idxs[end]
#         col_idxs = idxs[1:end-1]
#         col = LI[col_idxs...]
#         A[row, col] = val  # Assign the value directly
#     end

#     if nonredundant
#         return eliminate(A, p)  # Ensure 'eliminate' function is defined elsewhere
#     else
#         return A
#     end
# end

