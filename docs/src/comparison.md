# Product Comparison Guide

This page provides a comprehensive comparison of the different matrix products available in UniqueKronecker.jl.

## Overview of Products

UniqueKronecker.jl provides four main categories of products:

1. **Kronecker Product** (standard and unique)
2. **Khatri-Rao Product** (standard and unique)
3. **Face-Splitting Product** (standard and unique)
4. **Circulant Kronecker Product**

## Standard vs Unique Products

### Standard Products

Standard products include all possible products, including redundant terms:

```julia
using UniqueKronecker

x = [1, 2, 3]

# Standard Kronecker: includes both x[i]*x[j] and x[j]*x[i]
std_kron = kron(x, x)
# [1, 2, 3, 2, 4, 6, 3, 6, 9]
```

### Unique Products

Unique products eliminate redundancies by considering only unique index combinations:

```julia
# Unique Kronecker: includes only x[i]*x[j] where i ≤ j
unique_kron = x ⊘ x
# [1, 2, 3, 4, 6, 9]
```

## Product Comparison Table

| Product | Operator | Direction | Input Constraint | Output Size (for m×n input) |
|---------|----------|-----------|------------------|------------------------------|
| Kronecker | ⊗ | Element-wise | Vectors | ``n^2`` (standard), ``\frac{n(n+1)}{2}`` (unique) |
| Khatri-Rao | ⊙ | Column-wise | Same # columns | ``mp \times n`` (A: m×n, B: p×n) |
| Unique Khatri-Rao | ⨸ | Column-wise | Same # columns | Varies by unique terms per column |
| Face-Splitting | ⊖ | Row-wise | Same # rows | ``m \times np`` (A: m×n, B: m×p) |
| Unique Face-Splitting | ⧁ | Row-wise | Same # rows | Varies by unique terms per row |
| Circulant Kronecker | ⊛ | Cyclic permutations | Vectors/Matrices | Same as standard Kronecker ``n^k`` |

## Visual Examples

### Kronecker Product

```julia
x = [x₁, x₂]

# Standard
x ⊗ x = [x₁², x₁x₂, x₂x₁, x₂²]

# Unique
x ⊘ x = [x₁², x₁x₂, x₂²]
```

### Khatri-Rao Product (Column-wise)

```julia
A = [a₁ a₂]    B = [b₁ b₂]
    [a₃ a₄]        [b₃ b₄]

# Standard Khatri-Rao: A ⊙ B
Result = [a₁b₁  a₂b₂]
         [a₁b₃  a₂b₄]
         [a₃b₁  a₄b₂]
         [a₃b₃  a₄b₄]

# Column 1: kron([a₁,a₃], [b₁,b₃])
# Column 2: kron([a₂,a₄], [b₂,b₄])
```

### Face-Splitting Product (Row-wise)

```julia
A = [a₁ a₂]    B = [b₁ b₂]
    [a₃ a₄]        [b₃ b₄]

# Standard Face-Splitting: A ⊖ B
Result = [a₁b₁  a₁b₂  a₂b₁  a₂b₂]
         [a₃b₃  a₃b₄  a₄b₃  a₄b₄]

# Row 1: kron([a₁,a₂], [b₁,b₂])
# Row 2: kron([a₃,a₄], [b₃,b₄])
```

### Circulant Kronecker Product

The circulant Kronecker product sums over all cyclic permutations:

```julia
x = [x₁, x₂]
y = [y₁, y₂]

# Circulant Kronecker: x ⊛ y
x ⊛ y = (x ⊗ y) + (y ⊗ x)
      = [x₁y₁, x₁y₂, x₂y₁, x₂y₂] + [y₁x₁, y₁x₂, y₂x₁, y₂x₂]
      = [2x₁y₁, x₁y₂+y₁x₂, x₂y₁+y₂x₁, 2x₂y₂]

# For three vectors: x ⊛ y ⊛ z
x ⊛ y ⊛ z = (x ⊗ y ⊗ z) + (y ⊗ z ⊗ x) + (z ⊗ x ⊗ y)
```

## Use Case Selection Guide

### Choose Kronecker Product when:
- Working with single vectors
- Need element-wise products
- Building polynomial features for a single observation

### Choose Khatri-Rao Product when:
- Working with multiple feature vectors (columns)
- Each column represents a different sample/signal
- Tensor decompositions (CP/PARAFAC)
- Array signal processing

### Choose Face-Splitting Product when:

- Working with multiple observations (rows)
- Each row represents a different sample/time point
- Building polynomial features for multiple observations
- Dynamical systems with polynomial terms

### Choose Circulant Kronecker Product when:

- Computing derivatives of Kronecker products
- Working with symmetric polynomial structures
- Need all cyclic permutations of products
- Tensor calculus and differential operations
- Polynomial dynamical systems requiring symmetry

## Performance Considerations

### Memory Usage

For a vector ``\mathbf{x} \in \mathbb{R}^n``:

- Standard Kronecker (order ``k``): ``n^k`` elements
- Unique Kronecker (order ``k``): ``\binom{n+k-1}{k}`` elements

The reduction ratio for the unique product increases with dimension:

```julia
using UniqueKronecker

n = 10  # dimension
k = 2   # order

standard_size = n^k  # 100
unique_size = binomial(n + k - 1, k)  # 55

reduction = 1 - unique_size / standard_size  # 45% reduction
```

### Computational Efficiency

Unique products are generally faster when:
1. The dimension ``n`` is large
2. The order ``k`` is high
3. Subsequent computations only need unique terms

## Complete Example

```julia
using UniqueKronecker

# Example matrices
A = [1 2; 3 4]
B = [5 6; 7 8]
x = [2, 3, 4]

println("=== Kronecker Products ===")
kron_std = kron(x, x)
kron_uniq = x ⊘ x
println("Standard: ", kron_std)
println("Unique: ", kron_uniq)

println("\n=== Khatri-Rao Products (Column-wise) ===")
kr_std = A ⊙ B
kr_uniq = A ⨸ B
println("Standard size: ", size(kr_std))
println("Unique size: ", size(kr_uniq))

println("\n=== Face-Splitting Products (Row-wise) ===")
fs_std = A ⊖ B
fs_uniq = A ⧁ B
println("Standard size: ", size(fs_std))
println("Unique size: ", size(fs_uniq))

println("\n=== Powers ===")
# Khatri-Rao power
kr_power3 = A ⊙ 3
println("Khatri-Rao power 3 size: ", size(kr_power3))

# Face-splitting power
fs_power3 = A ⊖ 3
println("Face-splitting power 3 size: ", size(fs_power3))

println("\n=== Circulant Kronecker Products ===")
# Two vectors
y = [5, 6]
circ2 = x[1:2] ⊛ y
println("Circulant (2 vectors): ", circ2)

# Three vectors
z = [7, 8]
circ3 = x[1:2] ⊛ y ⊛ z
println("Circulant (3 vectors) size: ", size(circ3))
```

## Converting Between Standard and Unique

The package provides conversion functions:

```julia
using UniqueKronecker

# Create a standard Kronecker coefficient matrix
n = 3
A_std = zeros(n, n^2)
for i in 1:n
    x = rand(n)
    A_std[i,:] = kron(x, x)
end

# Convert to unique (eliminate redundancies)
A_uniq = eliminate(A_std, 2)

# Convert back to standard
A_back = duplicate(A_uniq, 2)

# Or with symmetric coefficients
A_sym = duplicate_symmetric(A_uniq, 2)
```

## Operator Summary

| Operator | LaTeX | How to Type | Function |
|----------|-------|-------------|----------|
| ⊗ | `\otimes` | `\\otimes<TAB>` | Standard Kronecker |
| ⊘ | `\oslash` | `\\oslash<TAB>` | Unique Kronecker |
| ⊙ | `\odot` | `\\odot<TAB>` | Khatri-Rao |
| ⨸ | - | `\\divideontimes<TAB>` | Unique Khatri-Rao |
| ⊖ | `\ominus` | `\\ominus<TAB>` | Face-Splitting |
| ⧁ | - | `\\revangle<TAB>` | Unique Face-Splitting |
| ⊛ | `\circledast` | `\\circledast<TAB>` | Circulant Kronecker |
