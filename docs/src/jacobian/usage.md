# Usage

## Quick Start

```julia
using UniqueKronecker

# Evaluate the Jacobian of the degree-2 unique monomial map at x = [3, 5]
x = [3.0, 5.0]
J = unique_kronecker_jacobian(x, 2)
# 3×2 Matrix{Float64}:
#  6.0  0.0
#  5.0  3.0
#  0.0  10.0
```

## Offline Precomputation with `DiffMatCache`

For repeated Jacobian evaluations (e.g., in a time-stepping loop or OpInf pipeline),
precompute the sparse differentiation matrices once and reuse them:

```julia
r = 5     # reduced dimension
pmax = 3  # maximum polynomial degree

# Precompute all N_d^{(i)} blocks for i = 1, ..., pmax
cache = DiffMatCache(r, pmax)

# Now evaluate Jacobians efficiently
x = randn(r)
J2 = unique_kronecker_jacobian(x, 2, cache)  # degree-2 Jacobian
J3 = unique_kronecker_jacobian(x, 3, cache)  # degree-3 Jacobian
```

## In-Place Evaluation

To avoid allocations in tight loops, use the in-place variant:

```julia
using LinearAlgebra: binomial

r = 4; deg = 3
cache = DiffMatCache(r, deg)

qi = binomial(r + deg - 1, deg)
J = zeros(qi, r)   # preallocate output

for x in eachcol(snapshot_matrix)
    unique_kronecker_jacobian!(J, x, deg, cache)
    # ... use J ...
end
```

## Batch Evaluation over Snapshot Matrices

Compute the Jacobian at every column of a snapshot matrix in one call:

```julia
r = 3; k = 100; deg = 2
X = randn(r, k)  # k snapshots of dimension r

# Returns a Vector{Matrix{Float64}} of length k
Js = unique_kronecker_jacobian(X, deg)

# With a precomputed cache
cache = DiffMatCache(r, deg)
Js = unique_kronecker_jacobian(X, deg, cache)
```

## Building Differentiation Matrices Directly

Access the sparse differentiation matrix blocks ``N_d^{(i)}`` for inspection
or custom computations:

```julia
# Individual blocks: N_1^{(2)}, N_2^{(2)} for r=2, i=2
Nblocks = diffmat_blocks(2, 2)
Nblocks[1]  # N_1^{(2)} — sparse 3×2
Nblocks[2]  # N_2^{(2)} — sparse 3×2

# Full concatenated matrix: N^{(i)} = [N_1^{(i)}, ..., N_r^{(i)}]
N = diffmat(2, 2)  # sparse 3×4
```

## Computing the Coupling Matrix

For the ManTA-OpInf augmented ROM, compute the coupling matrix
``C(\hat{s}) = \sum_{i=2}^{p} G^{(i)} \nabla_{\hat{s}} p_i^u(\hat{s})``:

```julia
r = 4; p = 3; j = 6

# Precompute differentiation matrices
cache = DiffMatCache(r, p)

# G^{(i)} = Vj' * V^{(i)} — precomputed projection matrices
G = [randn(j, binomial(r + i - 1, i)) for i in 2:p]

# Coupling matrix at a single point
x = randn(r)
C = coupling_matrix(x, G, cache)  # j × r matrix

# Batch: coupling matrix at each column of a snapshot matrix
X = randn(r, 100)
Cs = coupling_matrix(X, G, cache)  # Vector of j × r matrices

# The mass matrix for the augmented ROM
M = I(r) + C' * C   # symmetric positive definite
```

## Specialized Degree-2 and Degree-3 Jacobians

For maximum performance at low degrees, use the specialized routines that
bypass sparse matrix--vector products entirely:

```julia
x = randn(5)
J2 = unique_kronecker_jacobian2(x)  # degree-2, direct loop
J3 = unique_kronecker_jacobian3(x)  # degree-3, direct loop
```

These produce identical results to `unique_kronecker_jacobian(x, 2)` and
`unique_kronecker_jacobian(x, 3)`, respectively, but avoid sparse `mul!` overhead.

## Integration with the ManTA-OpInf Pipeline

Below is a sketch of how the Jacobian functions integrate into the full ManTA-OpInf 
Algorithm 1 (Formulation A):

```julia
using UniqueKronecker
using LinearAlgebra

# --- Phase 1–2: POD basis and polynomial manifold (assumed given) ---
# V   ∈ R^{n × r}        — POD basis
# Vbar ∈ R^{n × q}       — polynomial lifting operator [V^{(2)}, ..., V^{(p)}]
# Vj  ∈ R^{n × j}        — selected columns of Vbar
# Shat ∈ R^{r × k}       — reduced snapshots
# Sdothat ∈ R^{r × k}    — reduced time derivatives

# --- Phase 3: Offline precomputation for coupling matrix ---
r = size(V, 2)
cache = DiffMatCache(r, p)

# Precompute G^{(i)} = Vj' * V^{(i)}
G = Vector{Matrix{Float64}}(undef, p - 1)
col_offset = 0
for i in 2:p
    qi = binomial(r + i - 1, i)
    Vi = Vbar[:, col_offset+1 : col_offset+qi]
    G[i-1] = Vj' * Vi
    col_offset += qi
end

# Evaluate coupling matrices at each snapshot
Cs = coupling_matrix(Shat, G, cache)

# --- Phase 4: Operator inference (Formulation A) ---
# Assemble targets for the enriched block
Rbar = hcat([Cs[ℓ] * Sdothat[:, ℓ] for ℓ in 1:k]...)

# Solve the two decoupled OpInf problems using D, R, Rbar ...
```

---

## API Reference

### Differentiation Matrices

```@docs
diffmat_blocks
diffmat
```

### Cache

```@docs
DiffMatCache
```

### Jacobian Evaluation

```@docs
unique_kronecker_jacobian
unique_kronecker_jacobian!
unique_kronecker_jacobian2
unique_kronecker_jacobian3
```

### Coupling Matrix

```@docs
coupling_matrix
```