# Commutation Matrix

The commutation matrix $\mathbf{K}$ is a special permutation matrix that rearranges the vectorization of matrices and interchanges the Kronecker products of vectors and matrices. It plays a crucial role in simplifying expressions involving vectorization and Kronecker products in linear algebra, control theory, and other mathematical fields.

## Definition

The **Commutation Matrix** $\mathbf{K}_{m,n} \in \mathbb{R}^{mn \times mn}$ is defined as:

```math
\mathbf{K}_{m,n} = \sum_{i=1}^m \sum_{j=1}^n \mathbf{E}_{ij} \otimes \mathbf{E}_{ji}
```

where $\mathbf{E}_{ij}$ is the elementary matrix of size $m \times n$ with a one in the $(i, j)$-th position and zeros elsewhere.

**Properties of the Commutation Matrix:**

1. **Vectorization Transposition:**

   For any matrix $\mathbf{A} \in \mathbb{R}^{m \times n}$:

   ```math
   \mathbf{K}_{m,n} \operatorname{vec}(\mathbf{A}) = \operatorname{vec}(\mathbf{A}^\top)
   ```

   This property states that the commutation matrix transforms the vectorization of a matrix into the vectorization of its transpose.

2. **Kronecker Product Permutation:**

   For matrices $\mathbf{A} \in \mathbb{R}^{m \times n}$ and $\mathbf{B} \in \mathbb{R}^{p \times q}$:

   ```math
   \mathbf{K}_{p m, n q} (\mathbf{A} \otimes \mathbf{B}) \mathbf{K}_{n, q} = \mathbf{B} \otimes \mathbf{A}
   ```

   This property allows interchanging matrices within a Kronecker product.

3. **Vector Kronecker Product Permutation:**

   For vectors $\mathbf{v} \in \mathbb{R}^m$ and $\mathbf{w} \in \mathbb{R}^n$:

   ```math
   \mathbf{K}_{n,m} (\mathbf{v} \otimes \mathbf{w}) = \mathbf{w} \otimes \mathbf{v}
   ```

   This property enables swapping vectors within a Kronecker product.

**Additional Properties:**

- **Symmetry:** $\mathbf{K}_{n,m} = \mathbf{K}_{m,n}^\top$.
- **Involution:** $\mathbf{K}_{m,n} \mathbf{K}_{n,m} = \mathbf{I}_{mn}$, where $\mathbf{I}_{mn}$ is the identity matrix of size $mn \times mn$.
- **Eigenvalues:** The eigenvalues of $\mathbf{K}_{m,n}$ are either $+1$ or $-1$.

## Use Cases

The commutation matrix is particularly useful in the following contexts:

1. **Rearranging Kronecker Products:**

   The commutation matrix allows us to move matrices or vectors within Kronecker products. For example:

   ```math
   \mathbf{K}_{n,n^2} (\mathbf{A} \mathbf{x} \otimes \mathbf{x} \otimes \mathbf{x}) = \mathbf{x} \otimes \mathbf{A} \mathbf{x} \otimes \mathbf{x}
   ```

   ```math
   \mathbf{K}_{n^2,n} (\mathbf{A} \mathbf{x} \otimes \mathbf{x} \otimes \mathbf{x}) = \mathbf{x} \otimes \mathbf{x} \otimes \mathbf{A} \mathbf{x}
   ```

   Here, $\mathbf{x} \in \mathbb{R}^n$ and $\mathbf{A} \in \mathbb{R}^{n \times n}$.

2. **Simplifying Vectorized Equations:**

   In mathematical modeling, especially in statistics and econometrics, the commutation matrix simplifies expressions involving vectorized matrices.

3. **Matrix Derivatives:**

   It is used in deriving matrix derivatives where the order of differentiation and vectorization matters.

4. **Permutation of Tensor Products:**

   In higher-order tensor computations, the commutation matrix helps in permuting indices.

5. **Construct second-order symmetrizer matrix**

    For the second-order Kronecker product, the symmetrizer matrix is defined by

    ```math
        \mathbf{S}_2 = \frac{1}{2}(\mathbf{I}_{n^2} + \mathbf{K}_{n,n}) 
    ```

    See more details on the [symmetrizer matrix](symmetrizer.md).


```@docs
commat
```
