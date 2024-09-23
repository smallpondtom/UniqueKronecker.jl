# Unique Kronecker Powers

One unique feature of this package is the generalization of the unique Kronecker product for any *n*-dimensional vector with any Kronecker product of order *k*. This allows significant computational efficiency in higher-order polynomial systems.

## Definition

The **unique Kronecker power** of order $k$ of a vector $\mathbf{x} \in \mathbb{R}^n$, denoted by $\mathbf{x}^{\langle k \rangle}$, is a vector containing all the unique monomials of degree $k$ formed by the elements of $\mathbf{x}$. Unlike the standard Kronecker power $\mathbf{x}^{[k]}$, which includes all possible combinations (including duplicates due to commutativity), the unique Kronecker power includes each distinct monomial only once.

Formally, the unique Kronecker power is defined as:

```math
\mathbf{x}^{\langle k \rangle} = \mathrm{vec}_{\text{unique}} \left( \mathbf{x}^{[k]} \right) \in \mathbb{R}^{\tbinom{n + k - 1}{k}}
```

where:

- $\mathbf{x}^{[k]} = \underbrace{\mathbf{x} \otimes \mathbf{x} \otimes \cdots \otimes \mathbf{x}}_{k \text{ times}}$ is the $k$-th Kronecker power of $\mathbf{x}$.
- $\mathrm{vec}_{\text{unique}}$ extracts only the unique monomials, considering that the variables commute (e.g., $x_i x_j = x_j x_i$).

The length of $\mathbf{x}^{\langle k \rangle}$ is given by the multiset coefficient (number of combinations with repetition):

```math
\mathrm{length}(\mathbf{x}^{\langle k \rangle}) = \tbinom{n + k - 1}{k}
```

## Properties

- **Efficiency:** By considering only unique monomials, the unique Kronecker power significantly reduces the dimensionality compared to the standard Kronecker power, which has length $n^k$.
- **Symmetry:** The unique Kronecker power leverages the commutative property of multiplication, avoiding redundant terms.
- **Polynomial Representation:** Each element corresponds to a unique monomial of degree $k$ in the variables $x_1, x_2, \dots, x_n$.

## Computation

To compute the unique Kronecker power of a vector $\mathbf{x}$, generate all combinations of indices with replacement and compute the products corresponding to those combinations.

**Example for $n = 2$ and $k = 2$:**

Let $\mathbf{x} = [x_1, x_2]^\top$. The unique Kronecker power is:

```math
\mathbf{x}^{\langle 2 \rangle} = \begin{bmatrix} x_1^2 \\ x_1 x_2 \\ x_2^2 \end{bmatrix}
```

**Example for $n = 3$ and $k = 2$:**

```math
\mathbf{x}^{\langle 2 \rangle} = \begin{bmatrix} x_1^2 \\ x_1 x_2 \\ x_1 x_3 \\ x_2^2 \\ x_2 x_3 \\ x_3^2 \end{bmatrix}
```

## Conversion Between Standard and Unique Kronecker Powers

The package provides tools to convert between the standard Kronecker power $\mathbf{x}^{[k]}$ and the unique Kronecker power $\mathbf{x}^{\langle k \rangle}$ using the elimination matrix $\mathbf{L}_{n,k}$ and the duplication matrix $\mathbf{D}_{n,k}$.

**From Standard to Unique:**

```math
\mathbf{x}^{\langle k \rangle} = \mathbf{L}_{n,k} \mathbf{x}^{[k]}
```

**From Unique to Standard:**

```math
\mathbf{x}^{[k]} = \mathbf{D}_{n,k} \mathbf{x}^{\langle k \rangle}
```

```@docs
UniqueKronecker.UniqueCombinationIterator
```