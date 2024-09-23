# Duplication Matrix

In the second-order Kronecker product sense, the duplication matrix $\mathbf{D}$ is a matrix that transforms a vector containing only the unique elements of a symmetric matrix (vectorized in half-vectorization form) into the full vectorization of the symmetric matrix, effectively duplicating the symmetric elements to their appropriate positions. It plays a crucial role in manipulating Kronecker products and converting between non-redundant and redundant polynomial operators in polynomial dynamical systems.

## Definition

The duplication matrix $\mathbf{D}_n$ is a unique matrix of size $\left(n^2 \times \tfrac{n(n+1)}{2}\right)$ defined such that for any symmetric matrix $\mathbf{S} \in \mathbb{R}^{n \times n}$:

```math
    \mathrm{vec}(\mathbf{S}) = \mathbf{D}_n \mathrm{vech}(\mathbf{S})
```

where:

- $\mathrm{vec}(\mathbf{S})$ stacks all columns of $\mathbf{S}$ into a single vector of length $n^2$.
- $\mathrm{vech}(\mathbf{S})$ stacks only the elements on and below (or above) the main diagonal into a vector of length $\tfrac{n(n+1)}{2}$.

The duplication matrix effectively "duplicates" the unique elements in $\mathrm{vech}(\mathbf{S})$ to fill the positions in $\mathrm{vec}(\mathbf{S})$ that correspond to symmetric entries.

**Properties of the Duplication Matrix:**

- **Construction:** Each column of $\mathbf{D}_n$ contains ones at positions corresponding to the symmetric elements in $\operatorname{vec}(\mathbf{S})$ and zeros elsewhere.

- **Relation to Elimination Matrix:** The duplication matrix is related to the elimination matrix $\mathbf{L}_n$ via:

  ```math
      \mathbf{L}_n\mathbf{D}_n = \mathbf{I}_n
  ```

for more definitions see [MagnusEML1980](@citet).

!!! note
    Similar to the elimination matrix, this matrix is not formally/mathematically generalized for higher-order Kronecker products in the literature.

## Converting Between Kronecker Products

In polynomial dynamical systems, Kronecker powers of vectors often involve redundancy due to symmetric terms. The duplication matrix facilitates the conversion from the unique Kronecker product vector to the standard Kronecker product vector.

**Unique Kronecker Power:**

Given a vector $\mathbf{x} \in \mathbb{R}^n$, the unique Kronecker power is:

```math
    \mathbf{x}^{\langle k \rangle} \in \mathbb{R}^{\tbinom{n + k - 1}{k}}
```

which contains only the unique monomials of degree $k$.

**Standard Kronecker Power:**

The standard Kronecker power is:

```math
    \mathbf{x}^{[k]} = \underbrace{\mathbf{x} \otimes \mathbf{x} \otimes \cdots \otimes \mathbf{x}}_{k \text{ times}} \in \mathbb{R}^{n^k}
```

which includes all possible combinations, including redundant terms.

**Conversion Using the Duplication Matrix:**

The duplication matrix $\mathbf{D}_{n,k}$ of order $k$ is used to expand the unique Kronecker power to the standard Kronecker power:

```math
    \mathbf{x}^{[k]} = \mathbf{D}_{n,k} \mathbf{x}^{\langle k \rangle}
```

Here, $\mathbf{D}_{n,k}$ is a matrix of size $\left(n^k \times \tbinom{n + k - 1}{k}\right)$.

**Example:**

For $n = 2$ and $k = 2$, let $\mathbf{x} = [x_1, x_2]^\top$. Then:

- **Unique Kronecker Power:**

  ```math
      \mathbf{x}^{\langle 2 \rangle} = [x_1^2, x_1 x_2, x_2^2]^\top
  ```

- **Standard Kronecker Power:**

  ```math
      \mathbf{x}^{[2]} = [x_1^2, x_1 x_2, x_2 x_1, x_2^2]^\top
  ```

The duplication matrix $\mathbf{D}_{2,2}$ duplicates the mixed term $x_1 x_2$ to account for both $x_1 x_2$ and $x_2 x_1$ in the standard Kronecker power:

```math
    \mathbf{x}^{[2]} = \mathbf{D}_{2,2} \mathbf{x}^{\langle 2 \rangle}
```

**Note on Julia Implementation:**

The Julia function `dupmat(n::Int, p::Int)` generates the duplication matrix of order `p` for vectors of length `n`. This function can be used to create the duplication matrix needed for converting unique Kronecker products to standard ones.

```@docs
dupmat
```

**Benefits:**

- **Expands Unique Terms:** Allows reconstruction of the full Kronecker power from the unique monomials.
- **Simplifies Computations:** Facilitates operations that require the full Kronecker power.

## Converting Between Operators

The duplication matrix is also instrumental in converting between non-redundant (unique) and redundant polynomial operators in dynamical systems.

**Unique Operators:**

In the polynomial dynamical system:

```math
    \dot{\mathbf{x}}(t) = \mathbf{A} \mathbf{x}(t) + \mathbf{A}_{2u} \mathbf{x}^{\langle 2 \rangle}(t) + \mathbf{A}_{3u} \mathbf{x}^{\langle 3 \rangle}(t) + \cdots
```

the matrices $\mathbf{A}_{ku} \in \mathbb{R}^{n \times \tbinom{n + k - 1}{k}}$ operate on the unique Kronecker powers.

**Redundant Operators:**

Alternatively, using the standard Kronecker powers:

```math
    \dot{\mathbf{x}}(t) = \mathbf{A} \mathbf{x}(t) + \mathbf{A}_2 \mathbf{x}^{[2]}(t) + \mathbf{A}_3 \mathbf{x}^{[3]}(t) + \cdots
```

with $\mathbf{A}_k \in \mathbb{R}^{n \times n^k}$.

**Conversion Using the Duplication Matrix:**

The unique operator $\mathbf{A}_{ku}$ and the redundant operator $\mathbf{A}_k$ are related through the duplication matrix:

```math
    \mathbf{A}_{ku} = \mathbf{A}_{k} \mathbf{D}_{n,k}
```

which eliminates the redundant terms.


```@docs
eliminate
```