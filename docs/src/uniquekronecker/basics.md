# Unique Kronecker Product

## Introduction

The unique Kronecker product is a math operation which is similar to the Kronecker product but eliminates all the redundant terms appearing due to terms invariant under index permutations. This is easiest to understand through an example.

```@raw html
<img src="../../images/unique_kronecker_example.png" alt="unique kronecker example" style="float: left; margin-right: 10px;" />
```

The example above shows the unique Kronecker product for a vector of $\mathbf{x}\in\mathbb{R}^3$. For a standard Kronecker product the resulting vector becomes a 9-dimensional vector, and we have a 6-dimensional vector when using the unique Kronecker product. Which shows how the unique Kronecker product eliminates 3 of the redundant terms.

## Definitions

We define the basic syntax for the general Kronecker product and unique Kronecker product of order $k$ as follows:

```math
\begin{align}
    \mathbf{x}^{[k]} &= \underbrace{\mathbf{x} \otimes \cdots \otimes \mathbf{x}}_{k-\text{times}} \in \mathbb{R}^{n^k}, \\
    \mathbf{x}^{\langle k \rangle} &= \underbrace{\mathbf{x} \oslash \cdots \oslash \mathbf{x}}_{k-\text{times}} \in \mathbb{R}^{\binom{n+k-1}{k}}
\end{align}
```math

where $\mathbf{x}\in\mathbb{R}^n$. From this definition, we observe that the dimensions of the Kronecker product grows in the order of $n^k$ but the unique Kronecker grows in the order of $\binom{n+k-1}{k}$. The reduction in computational cost from this order difference becomes significantly obvious in higher-dimensions. 

For a second-order Kronecker product, the package supports the following syntax

```@repl
using UniqueKronecker
n = 3
x = rand(n)
x2u = x ⊘ x  # or unique_kronecker(x,x)
```

for higher-order Kronecker products, you can do

```@repl
using UniqueKronecker
n = 3
k = 5
x = rand(n)
xku = ⊘(x, k)
```

Another concept we need to define is the vectorization and half-vectorization operators $\mathrm{vec}(\cdot)$ and $\mathrm{vech}(\cdot)$, respectively. The vectorization operator flattens a matrix $\mathbf{A}\in\mathbb{R}^{m\times n}$ in the column direction creating a vector of size $mn$. On the other hand, the half-vectorization operator vectorizes the matrix but only half of it or discarding the supradiagonal entries. These operations are strongly related to the second-order Kronecker and unique Kronecker products, and those relationships are described in the picture below.

```@raw html
<img src="../../images/vectorize_and_kronecker.png" alt="unique kronecker example" style="float: left; margin-right: 10px;" />
```

The vectorization operator is already defined in the `LinearAlgebra` package as the function `vec`, but we define the functions `vech` for the half-vectorization operator and `invec` function for the inverse-vectorization operator that reverses the vectorization.

We are aware that similar concepts exists in the tensor algebra literature, and the vectorization and half-vectorization operations can be generalized to higher-order Kronecker products. However, for ease of exposition, we only illustrate the second-order Kronecker product case.


## Columnwise Operation on Snapshot Matrices

We also employ the functions 

- `kron_snapshot_matrix`
- `unique_kron_snapshot_matrix`

which allows you to apply the Kronecker product and unique Kronecker product on each column of a matrix. For example

```@repl
using UniqueKronecker
X = [1 2; 3 4]
X2 = kron_snapshot_matrix(X, 2)
```

```@repl
using UniqueKronecker
X = [1 2; 3 4]
X2u = unique_kron_snapshot_matrix(X, 2)
```


```@docs
unique_kronecker
⊘
kron_snapshot_matrix
unique_kron_snapshot_matrix
```
