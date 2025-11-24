# Khatri-Rao Product

## Overview

The **Khatri-Rao product** (denoted by ``\odot`` or `⊙` in this package) is a column-wise Kronecker product operation. Given two matrices ``\mathbf{A} \in \mathbb{R}^{m \times n}`` and ``\mathbf{B} \in \mathbb{R}^{p \times n}`` with the same number of columns, their Khatri-Rao product ``\mathbf{A} \odot \mathbf{B}`` is an ``mp \times n`` matrix where each column is the Kronecker product of the corresponding columns of ``\mathbf{A}`` and ``\mathbf{B}``:

```math
\mathbf{A} \odot \mathbf{B} = \begin{bmatrix} \mathbf{a}_1 \otimes \mathbf{b}_1 & \mathbf{a}_2 \otimes \mathbf{b}_2 & \cdots & \mathbf{a}_n \otimes \mathbf{b}_n \end{bmatrix}
```

where ``\mathbf{a}_i`` and ``\mathbf{b}_i`` are the ``i``-th columns of ``\mathbf{A}`` and ``\mathbf{B}``, respectively.

## Standard Khatri-Rao Product

### Basic Usage

```julia
using UniqueKronecker

A = [1 2; 3 4]
B = [5 6; 7 8]

# Standard Khatri-Rao product
kr = khatri_rao(A, B)
# Equivalent to:
# Column 1: kron([1,3], [5,7]) = [5, 7, 15, 21]
# Column 2: kron([2,4], [6,8]) = [12, 16, 24, 32]
```

### Unicode Operator

You can use the Unicode operator `⊙`:

```julia
kr = A ⊙ B  # Same as khatri_rao(A, B)
```

### Self-Product

When called with a single matrix argument, it computes the Khatri-Rao product of the matrix with itself:

```julia
kr_self = khatri_rao(A)  # Same as khatri_rao(A, A)
kr_self = A ⊙ A          # or using the operator
```

### Multiple Matrices

The Khatri-Rao product can be extended to multiple matrices:

```julia
A = [1 2; 3 4]
B = [5 6; 7 8]
C = [9 10; 11 12]

kr_multi = khatri_rao(A, B, C)
kr_multi = A ⊙ B ⊙ C  # Using operators
```

### Power Form

You can compute the ``d``-th Khatri-Rao power of a matrix, where each column is raised to the Kronecker power ``d``:

```julia
A = [1 2; 3 4]

# Second Khatri-Rao power
kr_power2 = khatri_rao(A, 2)  # Each column: kron(A[:,j], A[:,j])
kr_power2 = A ⊙ 2

# Third Khatri-Rao power
kr_power3 = khatri_rao(A, 3)  # Each column: kron(kron(A[:,j], A[:,j]), A[:,j])
kr_power3 = A ⊙ 3
```

## Unique Khatri-Rao Product

The **unique Khatri-Rao product** (denoted by `⨸` in this package) eliminates redundant terms in each column by applying the unique Kronecker product to each column pair:

```math
\mathbf{A}\,⨸\,\mathbf{B} = \begin{bmatrix} \mathbf{a}_1 \oslash \mathbf{b}_1 & \mathbf{a}_2 \oslash \mathbf{b}_2 & \cdots & \mathbf{a}_n \oslash \mathbf{b}_n \end{bmatrix}
```

### Basic Usage

```julia
using UniqueKronecker

A = [1 2; 3 4]
B = [5 6; 7 8]

# Unique Khatri-Rao product
ukr = unique_khatri_rao(A, B)
ukr = A ⨸ B  # Using the operator

# For self-product with redundancies eliminated
ukr_self = unique_khatri_rao(A)
ukr_self = A ⨸ A
```

### Power Form

```julia
A = [1 2; 3 4]

# Second unique Khatri-Rao power
ukr_power2 = unique_khatri_rao(A, 2)
ukr_power2 = A ⨸ 2

# Third unique Khatri-Rao power
ukr_power3 = unique_khatri_rao(A, 3)
ukr_power3 = A ⨸ 3
```

## Applications

### Tensor Decompositions

The Khatri-Rao product is fundamental in tensor decompositions such as CANDECOMP/PARAFAC (CP) decomposition:

```julia
# In CP decomposition, the tensor approximation involves Khatri-Rao products
# of factor matrices
A = rand(10, 3)  # Factor matrix 1
B = rand(15, 3)  # Factor matrix 2
C = rand(20, 3)  # Factor matrix 3

# The unfolded tensor approximation
approx = (C ⊙ B ⊙ A)
```

### Signal Processing

In array signal processing, the Khatri-Rao product appears in the modeling of sensor array outputs:

```julia
# Steering vectors for different arrays
A = rand(8, 4)   # Array 1 steering vectors
B = rand(6, 4)   # Array 2 steering vectors

# Combined array response
response = A ⊙ B
```

### Machine Learning Feature Engineering

For polynomial feature expansion with structured interactions:

```julia
# Feature matrices from different modalities
X1 = rand(100, 5)  # Features from modality 1
X2 = rand(100, 5)  # Features from modality 2

# Create interaction features column-wise
interactions = unique_khatri_rao(X1, X2)
```

## Properties

### Size

For ``\mathbf{A} \in \mathbb{R}^{m \times n}`` and ``\mathbf{B} \in \mathbb{R}^{p \times n}``:
- Standard Khatri-Rao: ``\mathbf{A} \odot \mathbf{B} \in \mathbb{R}^{mp \times n}``
- Unique Khatri-Rao (when ``m = p``): Result has ``\frac{m(m+1)}{2}`` rows for the self-product

### Relationship to Kronecker Product

The Khatri-Rao product can be viewed as a column-wise application of the Kronecker product:

```julia
A = [1 2; 3 4]
B = [5 6; 7 8]

# Khatri-Rao
kr = A ⊙ B

# Manual construction
kr_manual = hcat([kron(A[:,i], B[:,i]) for i in 1:size(A,2)]...)

@assert kr == kr_manual
```

## API Reference

```@docs
khatri_rao
unique_khatri_rao
⊙
⨸
```
