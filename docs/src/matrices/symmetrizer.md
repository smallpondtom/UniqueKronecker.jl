# Symmetrizer Matrix

The symmetrizer matrix $\mathbf{S}$ is a matrix that symmetrizes tensors or higher-order Kronecker products. It ensures that the resulting product is symmetric with respect to its indices. In the context of polynomial dynamical systems and Kronecker products, the symmetrizer matrix is crucial for handling symmetric terms efficiently and accurately.

## Definition

The **Symmetrizer Matrix** $\mathbf{S}_{n,p} \in \mathbb{R}^{n^p \times n^p}$ of order $p$ for a vector of length $n$ is defined such that it symmetrizes a tensor obtained from the $p$-th Kronecker power of a vector $\mathbf{x} \in \mathbb{R}^n$.

Given the $p$-th Kronecker power:

```math
\mathbf{x}^{[p]} = \underbrace{\mathbf{x} \otimes \mathbf{x} \otimes \cdots \otimes \mathbf{x}}_{p \text{ times}} \in \mathbb{R}^{n^p}
```

The symmetrizer matrix $\mathbf{S}_{n,p}$ symmetrizes $\mathbf{x}^{[p]}$ to produce a vector where all permutations of indices are averaged. The symmetrized vector is:

```math
\mathbf{x}^{\{p\}} = \mathbf{S}_{n,p} \mathbf{x}^{[p]}
```

**Properties of the Symmetrizer Matrix:**

- **Idempotent:** $\mathbf{S}_{n,p}^2 = \mathbf{S}_{n,p}$.
- **Symmetric:** $\mathbf{S}_{n,p}^\top = \mathbf{S}_{n,p}$.
- **Projection Matrix:** It projects any tensor onto the space of symmetric tensors.

For the quadratic case ($p = 2$), the symmetrizer matrix can be explicitly defined using the commutation matrix $\mathbf{K}_{n,n}$:

```math
\mathbf{S}_{n,2} = \frac{1}{2} \left( \mathbf{I}_{n^2} + \mathbf{K}_{n,n} \right)
```

where $\mathbf{I}_{n^2}$ is the identity matrix of size $n^2 \times n^2$.

```@docs
symmtzrmat
```

## Use Cases

The symmetrizer matrix is particularly useful in the following contexts:

1. **Symmetrizing Kronecker Products:**

   When dealing with polynomial terms in dynamical systems, it's often necessary to ensure that the terms are symmetric. The symmetrizer matrix averages all permutations of indices in the Kronecker product to produce a symmetric tensor.

   **Example:**

   For $\mathbf{x} \in \mathbb{R}^n$ and $p = 2$:

   ```math
   \mathbf{x}^{\{2\}} = \mathbf{S}_{n,2} (\mathbf{x} \otimes \mathbf{x})
   ```

   This results in a vector containing terms like $x_i x_j$ averaged over all permutations.

2. **Ensuring Symmetric Operators:**

   In polynomial dynamical systems, operators (matrices) that act on Kronecker powers of state vectors should often be symmetric to preserve certain properties like stability. The symmetrizer matrix can be used to enforce symmetry in these operators.

3. **Duplicating Symmetric Coefficients:**

   When converting between unique and redundant representations of polynomial operators, the symmetrizer matrix ensures that coefficients are duplicated symmetrically, maintaining the symmetry of the operator.

```docs
duplicate_symmetric
```