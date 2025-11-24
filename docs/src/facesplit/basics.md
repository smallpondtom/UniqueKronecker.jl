# Face-Splitting Product

## Overview

The **face-splitting product** (denoted by `⊖` in this package) is a row-wise Kronecker product operation. Given two matrices ``\mathbf{A} \in \mathbb{R}^{m \times n}`` and ``\mathbf{B} \in \mathbb{R}^{m \times p}`` with the same number of rows, their face-splitting product ``\mathbf{A} \ominus \mathbf{B}`` is an ``m \times np`` matrix where each row is the Kronecker product of the corresponding rows of ``\mathbf{A}`` and ``\mathbf{B}``:

```math
\mathbf{A} \ominus \mathbf{B} = \begin{bmatrix} 
(\mathbf{a}_1^\top \otimes \mathbf{b}_1^\top)^\top \\[0.2em] 
(\mathbf{a}_2^\top \otimes \mathbf{b}_2^\top)^\top \\[0.2em]
\vdots \\[0.2em]
(\mathbf{a}_m^\top \otimes \mathbf{b}_m^\top)^\top 
\end{bmatrix}
```

where ``\mathbf{a}_i^\top`` and ``\mathbf{b}_i^\top`` are the ``i``-th rows of ``\mathbf{A}`` and ``\mathbf{B}``, respectively.

## Standard Face-Splitting Product

### Usage Examples

```julia
using UniqueKronecker

A = [1 2; 3 4]
B = [5 6; 7 8]

# Standard face-splitting product
fs = face_split(A, B)
# Equivalent to:
# Row 1: kron([1,2], [5,6]) = [5, 6, 10, 12]
# Row 2: kron([3,4], [7,8]) = [21, 24, 28, 32]
```

### Unicode Operator

You can use the Unicode operator `⊖`:

```julia
fs = A ⊖ B  # Same as face_split(A, B)
```

### Self-Product

When called with a single matrix argument, it computes the face-splitting product of the matrix with itself:

```julia
fs_self = face_split(A)  # Same as face_split(A, A)
fs_self = A ⊖ A          # or using the operator
```

### Multiple Matrices

The face-splitting product can be extended to multiple matrices:

```julia
A = [1 2; 3 4]
B = [5 6; 7 8]
C = [9 10; 11 12]

fs_multi = face_split(A, B, C)
fs_multi = A ⊖ B ⊖ C  # Using operators
```

### Power Form

You can compute the ``d``-th face-splitting power of a matrix, where each row is raised to the Kronecker power ``d``:

```julia
A = [1 2; 3 4]

# Second face-splitting power
fs_power2 = face_split(A, 2)  # Each row: kron(A[i,:], A[i,:])
fs_power2 = A ⊖ 2

# Third face-splitting power
fs_power3 = face_split(A, 3)  # Each row: kron(kron(A[i,:], A[i,:]), A[i,:])
fs_power3 = A ⊖ 3
```

## Unique Face-Splitting Product

The **unique face-splitting product** (denoted by `⧁` in this package) eliminates redundant terms in each row by applying the unique Kronecker product to each row pair:

```math
\mathbf{A}\,⧁\,\mathbf{B} = \begin{bmatrix} 
(\mathbf{a}_1^\top \oslash \mathbf{b}_1^\top)^\top \\[0.2em]
(\mathbf{a}_2^\top \oslash \mathbf{b}_2^\top)^\top \\[0.2em]
\vdots \\[0.2em]
(\mathbf{a}_m^\top \oslash \mathbf{b}_m^\top)^\top 
\end{bmatrix}
```

### Usage Examples

```julia
using UniqueKronecker

A = [1 2; 3 4]
B = [5 6; 7 8]

# Unique face-splitting product
ufs = unique_face_split(A, B)
ufs = A ⧁ B  # Using the operator

# For self-product with redundancies eliminated
ufs_self = unique_face_split(A)
ufs_self = A ⧁ A
```

### Power Form with Unique Products

```julia
A = [1 2; 3 4]

# Second unique face-splitting power
ufs_power2 = unique_face_split(A, 2)
ufs_power2 = A ⧁ 2

# Third unique face-splitting power
ufs_power3 = unique_face_split(A, 3)
ufs_power3 = A ⧁ 3
```

## Applications

### Matrix Polynomial Operations

The face-splitting product is useful for representing polynomial operations on matrix rows:

```julia
# Data matrix
X = rand(100, 5)

# Create quadratic features row-wise (each sample independently)
X_quad = unique_face_split(X, 2)

# Each row of X_quad contains unique quadratic terms for that sample
```

### Dynamical Systems

In polynomial dynamical systems, face-splitting products can represent state transformations:

```julia
# State trajectory matrix (each row is a time point)
X = rand(50, 3)

# Create polynomial features for each time point
X_poly = unique_face_split(X, 2)

# Use in dynamical system identification
# X_dot = X_poly * coefficients
```

### Structured Matrix Computations

For operations requiring row-wise polynomial expansions:

```julia
# Input data with multiple observations
A = rand(10, 4)

# Create second-order row-wise features
A2 = face_split(A, 2)

# Each row of A2 contains all pairwise products for that observation
```

## Properties

### Size

For ``\mathbf{A} \in \mathbb{R}^{m \times n}`` and ``\mathbf{B} \in \mathbb{R}^{m \times p}``:

- Standard face-splitting: ``\mathbf{A} \ominus \mathbf{B} \in \mathbb{R}^{m \times np}``
- Unique face-splitting (when ``n = p``): Result has ``\frac{n(n+1)}{2}`` columns for the self-product

### Relationship to Kronecker Product

The face-splitting product can be viewed as a row-wise application of the Kronecker product:

```julia
A = [1 2; 3 4]
B = [5 6; 7 8]

# Face-splitting
fs = A ⊖ B

# Manual construction
fs_manual = vcat([transpose(kron(A[i,:], B[i,:])) for i in 1:size(A,1)]...)

@assert fs == fs_manual
```

## Comparison: Khatri-Rao vs Face-Splitting

| Property | Khatri-Rao (⊙) | Face-Splitting (⊖) |
|----------|----------------|-------------------|
| Operation | Column-wise | Row-wise |
| Input constraint | Same # of columns | Same # of rows |
| Output for ``m \times n`` ⊙/⊖ ``p \times q`` | ``mp \times n`` (if ``n = q``) | ``m \times np`` (if ``m = p``) |
| Primary use | Tensor decompositions | Polynomial features |
| Each element processes | Columns (features) | Rows (observations) |

## API Reference

```@docs
face_split
unique_face_split
⊖
⧁
```
