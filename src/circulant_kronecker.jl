export circulant_kronecker, ⊛, circulant_kron_snapshot_matrix


"""
    circulant_kronecker(args::AbstractArray...)

Circulant Kronecker product operation for multiple Kronecker products.

## Arguments
- `args::AbstractArray...`: Vectors or matrices to perform the circulant Kronecker product.

## Returns
- `result`: The circulant Kronecker product.
"""
function circulant_kronecker(args::AbstractArray...)
    n = length(args)
    if n == 2
        return args[1] ⊗ args[2] + args[2] ⊗ args[1]
    elseif n == 3
        return args[1] ⊗ args[2] ⊗ args[3] + args[2] ⊗ args[3] ⊗ args[1] + args[3] ⊗ args[1] ⊗ args[2]
    else
        # Compute one Kronecker product to get the correct size and type
        initial_kp = ⊗(args...)
        # Initialize result as a zero array of the same type and size as initial_kp
        result = zero(initial_kp)
        @inbounds for k in 0:(n - 1)
            # Generate cyclic permutation of arguments
            permuted_args = ntuple(i -> args[mod1(i + k, n)], n)
            # Compute Kronecker product of permuted arguments
            kp = ⊗(permuted_args...)
            result += kp
        end
        return result
    end
end

# Helper function to flatten arguments
function flatten_args(args...)
    result = []
    for arg in args
        if isa(arg, AbstractArray{<:Number})
            push!(result, arg)
        elseif isa(arg, Tuple) || isa(arg, Array)
            append!(result, flatten_args(arg...))
        else
            error("Unsupported argument type: ", typeof(arg))
        end
    end
    return result
end

"""
    ⊛(args::AbstractArray...)

Circulant Kronecker product operator for multiple arguments.

## Arguments
- `args::AbstractArray...`: Vectors or matrices to perform the circulant Kronecker product.

## Returns
- `result`: The circulant Kronecker product.
"""
function ⊛(args...)
    flat_args = flatten_args(args...)
    return circulant_kronecker(flat_args...)
end


"""
    circulant_kron_snapshot_matrix(Xmat::AbstractArray{T}...) where {T<:Number}

Compute the circulant Kronecker product of a set of matrices, where each matrix is a snapshot matrix.

## Arguments
- `Xmat::AbstractArray{T}...`: Snapshot matrices to compute the circulant Kronecker product.

## Returns
- `result`: The circulant Kronecker product of the snapshot matrices.
"""
function circulant_kron_snapshot_matrix(Xmat::AbstractArray{T}...) where {T<:Number}
    # Ensure all matrices have the same number of columns
    ncols = size(Xmat[1], 2)
    for X in Xmat[2:end]
        if size(X, 2) != ncols
            throw(ArgumentError("All input matrices must have the same number of columns"))
        end
    end

    # Inner function to compute circulant Kronecker product for a set of vectors
    function circulant_kron_timestep(x...)
        return circulant_kronecker(x...)
    end

    # Collect columns from each matrix
    col_iterators = map(eachcol, Xmat)
    # Zip the columns together to align them
    zipped_columns = zip(col_iterators...)
    # Apply circulant_kron_timestep to each set of columns
    tmp = [circulant_kron_timestep(cols...) for cols in zipped_columns]
    # Concatenate the results horizontally
    return reduce(hcat, tmp)
end
