# Theory

## Jacobian of the Unique Monomial Map

Given a vector ``\hat{s} \in \mathbb{R}^r``, the unique monomial map of degree ``i`` is the vector

```math
p_i^u(\hat{s}) \in \mathbb{R}^{q_i}, \quad q_i := \binom{r + i - 1}{i},
```

whose entries are all unique monomials of total degree ``i`` in the components of ``\hat{s}``,
ordered in non-decreasing index order.  For example, with ``r = 2``:

```math
p_2^u(\hat{s}) = \begin{bmatrix} \hat{s}_1^2 \\ \hat{s}_1 \hat{s}_2 \\ \hat{s}_2^2 \end{bmatrix}, \qquad
p_3^u(\hat{s}) = \begin{bmatrix} \hat{s}_1^3 \\ \hat{s}_1^2 \hat{s}_2 \\ \hat{s}_1 \hat{s}_2^2 \\ \hat{s}_2^3 \end{bmatrix}.
```

The Jacobian ``\nabla_{\hat{s}} p_i^u(\hat{s}) \in \mathbb{R}^{q_i \times r}`` is the matrix of partial
derivatives of these monomials with respect to ``\hat{s}``.

## Differentiation Matrices

Rather than computing the Jacobian via the full Kronecker product, we use a factored form involving
sparse **differentiation matrices** ``N_d^{(i)} \in \mathbb{R}^{q_i \times q_{i-1}}`` for
``d = 1, \dots, r``.

Each differentiation matrix is defined entry-wise by

```math
[N_d^{(i)}]_{\alpha,\beta} =
\begin{cases}
\alpha_d & \text{if } \beta = \alpha - e_d \text{ and } \alpha_d \ge 1, \\
0 & \text{otherwise},
\end{cases}
```

where ``\alpha \in \mathcal{I}_i`` and ``\beta \in \mathcal{I}_{i-1}`` are multi-indices with ``|\alpha| = i``
and ``|\beta| = i-1``, and ``e_d`` is the ``d``-th coordinate unit multi-index.

!!! info "Key properties"
    - Each ``N_d^{(i)}`` has **at most one nonzero per row**, making it extremely sparse.
    - The total number of nonzeros across all ``r`` blocks is ``r \cdot \binom{r+i-2}{i-1}``.
    - The sum of all values across all blocks equals ``i \cdot q_i``.
    - All matrices are **constant** (independent of ``\hat{s}``), so they can be precomputed once and stored.

The Jacobian formula then reads:

```math
\nabla_{\hat{s}} p_i^u(\hat{s}) = N^{(i)} \bigl(I_r \otimes p_{i-1}^u(\hat{s})\bigr),
```

where ``N^{(i)} = [N_1^{(i)}, \dots, N_r^{(i)}] \in \mathbb{R}^{q_i \times r \, q_{i-1}}``
is the horizontal concatenation. Equivalently, the ``d``-th column of the Jacobian is

```math
\frac{\partial p_i^u(\hat{s})}{\partial \hat{s}_d} = N_d^{(i)} \, p_{i-1}^u(\hat{s}),
```

which is simply a sparse matrix--vector product.

## Worked Example (``r = 2``, ``i = 2``)

The index sets are ``\mathcal{I}_2 = \{(2,0), (1,1), (0,2)\}`` (rows) and
``\mathcal{I}_1 = \{(1,0), (0,1)\}`` (columns).

**Differentiation matrix** ``N_1^{(2)}``:

| Row ``\alpha`` | ``\alpha - e_1`` | ``\alpha_1`` | Column |
|:---:|:---:|:---:|:---:|
| ``(2,0)`` | ``(1,0)`` | ``2`` | ``1`` |
| ``(1,1)`` | ``(0,1)`` | ``1`` | ``2`` |
| ``(0,2)`` | — | ``0`` | — |

```math
N_1^{(2)} = \begin{bmatrix} 2 & 0 \\ 0 & 1 \\ 0 & 0 \end{bmatrix}
```

**Differentiation matrix** ``N_2^{(2)}``:

| Row ``\alpha`` | ``\alpha - e_2`` | ``\alpha_2`` | Column |
|:---:|:---:|:---:|:---:|
| ``(2,0)`` | — | ``0`` | — |
| ``(1,1)`` | ``(1,0)`` | ``1`` | ``1`` |
| ``(0,2)`` | ``(0,1)`` | ``2`` | ``2`` |

```math
N_2^{(2)} = \begin{bmatrix} 0 & 0 \\ 1 & 0 \\ 0 & 2 \end{bmatrix}
```

**Jacobian evaluation** at ``\hat{s} = [3, 5]^\top``:

```math
\nabla_{\hat{s}} p_2^u(\hat{s}) =
\begin{bmatrix}
N_1^{(2)} p_1^u(\hat{s}) & N_2^{(2)} p_1^u(\hat{s})
\end{bmatrix} =
\begin{bmatrix}
2 \cdot 3 & 0 \\
5 & 3 \\
0 & 2 \cdot 5
\end{bmatrix} =
\begin{bmatrix}
6 & 0 \\
5 & 3 \\
0 & 10
\end{bmatrix}
```

which matches the expected derivative of ``p_2^u(\hat{s}) = [\hat{s}_1^2, \hat{s}_1\hat{s}_2, \hat{s}_2^2]^\top``.

## Coupling Matrix

In the ManTA-OpInf framework, the **coupling matrix** 
``C(\hat{s}) \in \mathbb{R}^{j \times r}`` links the standard and enriched 
operator inference problems:

```math
C(\hat{s}) = \sum_{i=2}^{p} G^{(i)} \nabla_{\hat{s}} p_i^u(\hat{s}),
```

where ``G^{(i)} = V_j^\top V^{(i)} \in \mathbb{R}^{j \times q_i}`` are precomputed projection
matrices. The coupling matrix appears in both the mass matrix 
``M(\hat{s}) = I_r + C(\hat{s})^\top C(\hat{s})`` and the augmented right-hand side of the 
reduced-order model.

## Relation to the Elimination Matrix

The differentiation matrices provide an alternative to the route via the elimination and symmetrizer matrix 
``D_{r,i}``. While the identity ``\nabla_{\hat{s}} p_i^u(\hat{s}) = L_{r,i} S_{r,i} \nabla_{\hat{s}} (\hat{s}^{\otimes i})``
is mathematically equivalent, the approach via ``N_d^{(i)}`` avoids forming the full 
``r^i``-dimensional Kronecker Jacobian entirely and exploits the fact that 
``p_{i-1}^u(\hat{s})`` is already computed during ROM feature evaluation.