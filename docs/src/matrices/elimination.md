# Elimination Matrix

Elimination matrix $\mathbf{L}$, first mentioned by Tracy and Singh (1972) and later by Vetter (1975), was studied in depth by Magnus and Neudecker in [MagnusEML1980](@citet). This matrix holds properties that are crucial to converting the standard Kronecker product vector to the unique Kronecker product vector as well as converting the non-redundant polynomial operators to the redundant operators.

## Definition

The elimination matrix $\mathbf{L}_n$ is a unique matrix of size $\left(\tfrac{n(n+1)}{2} \times n^2\right)$ that extracts the unique elements from the vectorization of a symmetric matrix. For any symmetric matrix $\mathbf{S} \in \mathbb{R}^{n \times n}$, the vectorization $\mathrm{vec}(\mathbf{S})$ contains redundant elements due to symmetry. The elimination matrix $\mathbf{L}_n$ is defined such that:

```math
    \mathrm{vech}(\mathbf{S}) = \mathbf{L}_n \mathrm{vec}(\mathbf{S})
```

where:

- operator $\mathrm{vec}(\mathbf{S})$ stacks all columns of $\mathbf{S}$ into a single vector of length $n^2$.
- and $\mathrm{vech}(\mathbf{S})$ stacks only the elements on and below (or above) the main diagonal into a vector of length $\tfrac{n(n+1)}{2}$.

The elimination matrix effectively "eliminates" the redundant elements in $\mathrm{vec}(\mathbf{S})$, preserving only the unique components.

**Properties of the Elimination Matrix:**

- **Idempotency:** $\mathbf{L}_n \mathbf{L}_n^\top$ is an idempotent matrix, meaning $(\mathbf{L}_n \mathbf{L}_n^\top)^2 = \mathbf{L}_n \mathbf{L}_n^\top$.
- **Construction:** Each row of $\mathbf{L}_n$ contains a single one and zeros elsewhere, mapping elements from $\mathrm{vec}(\mathbf{S})$ to $\mathrm{vech}(\mathbf{S})$.

!!! note
    Unfortunately, the definition of the elimination matrix is not formally/mathematically defined for Kronecker products beyond the order of 2. However, implementation-wise, this package has generalized the matrix.

## Converting Between Kronecker Products

In polynomial dynamical systems, we often deal with Kronecker powers of vectors, which can be redundant due to symmetric terms. The elimination matrix facilitates the conversion between the standard Kronecker product vector and the unique Kronecker product vector.

**Standard Kronecker Power:**

Given a vector $\mathbf{x} \in \mathbb{R}^n$, the $k$-th order Kronecker power is:

```math
    \mathbf{x}^{[k]} = \underbrace{\mathbf{x} \otimes \mathbf{x} \otimes \cdots \otimes \mathbf{x}}_{k \text{ times}} \in \mathbb{R}^{n^k}
```

**Unique Kronecker Power:**

The unique Kronecker power, denoted by $\mathbf{x}^{\langle k \rangle}$, contains only the unique monomials of degree $k$ formed from the elements of $\mathbf{x}$. It has a length of $\tbinom{n + k - 1}{k}$.

**Conversion Using the Elimination Matrix:**

The relationship between the standard and unique Kronecker powers is given by:

```math
    \mathbf{x}^{\langle k \rangle} = \mathbf{L}_{n,k} \mathbf{x}^{[k]}
```

where $\mathbf{L}_{n,k}$ is the elimination matrix corresponding to the $k$-th Kronecker power, of size $\left(\tbinom{n + k - 1}{k} \times n^k\right)$.

**Example:**

For $n = 2$ and $k = 2$, let $\mathbf{x} = [x_1, x_2]^\top$. Then:

- **Standard Kronecker Power:**

  ```math
      \mathbf{x}^{[2]} = [x_1^2, x_1 x_2, x_2 x_1, x_2^2]^\top
  ```

- **Unique Kronecker Power:**

  ```math
      \mathbf{x}^{\langle 2 \rangle} = [x_1^2, x_1 x_2, x_2^2]^\top
  ```

The elimination matrix $\mathbf{L}_{2,2}$ maps $\mathbf{x}^{[2]}$ to $\mathbf{x}^{\langle 2 \rangle}$ by combining redundant terms (since $x_1 x_2$ and $x_2 x_1$ are the same):

```math
    \mathbf{x}^{\langle 2 \rangle} = \mathbf{L}_{2,2} \mathbf{x}^{[2]}
```


```@docs
elimat
```

**Benefits:**

- **Efficiency:** Reduces computational complexity by working with smaller vectors.
- **Simplicity:** Simplifies expressions and computations involving symmetric terms.
- **Clarity:** Makes the structure of polynomial terms more transparent.

## Converting Between Operators

The elimination matrix also plays a crucial role in converting between redundant and non-redundant polynomial operators in dynamical systems.

**Redundant Operators:**

In the polynomial system:

```math
    \dot{\mathbf{x}}(t) = \mathbf{A}\mathbf{x}(t) + \mathbf{A}_2 \mathbf{x}^{[2]}(t) + \mathbf{A}_3 \mathbf{x}^{[3]}(t) + \cdots
```

the matrices $\mathbf{A}_k \in \mathbb{R}^{n \times n^k}$ operate on the full Kronecker powers, which include redundant terms.

**Unique Operators:**

To eliminate redundancy, we define unique operators $\mathbf{A}_{ku} \in \mathbb{R}^{n \times \tbinom{n + k - 1}{k}}$ such that:

```math
    \dot{\mathbf{x}}(t) = \mathbf{A}\mathbf{x}(t) + \mathbf{A}_{2u} \mathbf{x}^{\langle 2 \rangle}(t) + \mathbf{A}_{3u} \mathbf{x}^{\langle 3 \rangle}(t) + \cdots
```

**Conversion Using the Elimination Matrix:**

The redundant operator $\mathbf{A}_k$ and the unique operator $\mathbf{A}_{ku}$ are related through the elimination matrix:

```math
    \mathbf{A}_{k} = \mathbf{A}_{ku} \mathbf{L}_{n,k}
```

which actually works to duplicate the entries of the non-redundant operator. Thus, though seemingly counter-intuitive, the elimination matrix duplicates the entries and constructs to convert the $\mathbf{A}_2$ matrix to the $\mathbf{A}_{2u}$ matrix.


**Example:**

Consider a system with $n = 2$:

- **Redundant operator:** $\mathbf{A}_2 \in \mathbb{R}^{2 \times 4}$
- **Unique operator:** $\mathbf{A}_{2u} \in \mathbb{R}^{2 \times 3}$

The elimination matrix $\mathbf{L}_{2,2}$ maps $\mathbf{x}^{[2]}$ to $\mathbf{x}^{\langle 2 \rangle}$:

```math
    \mathbf{x}^{\langle 2 \rangle} = \mathbf{L}_{2,2} \mathbf{x}^{[2]}
```

The operators are related by:

```math
    \mathbf{A}_{2u} = \mathbf{A}_2 \mathbf{L}_{2,2}^\top
```

```@docs
duplicate
```