# Circulant Kronecker Product

The **Circulant Kronecker Product** is a generalized operation that extends the standard Kronecker product by summing over all cyclic permutations of the given vectors or matrices. This operation can be applied to any number of vectors or matrices.

## Mathematical Definition

Given a sequence of vectors or matrices $\mathbf{x}_1, \mathbf{x}_2, \dots, \mathbf{x}_n$, the circulant Kronecker product is defined as:

```math
\mathbf{x}_1 \circledast \mathbf{x}_2 \circledast \dots \circledast \mathbf{x}_n = \sum_{k=0}^{n-1} \mathbf{x}_{1+k} \otimes \mathbf{x}_{2+k} \otimes \dots \otimes \mathbf{x}_{n+k}
```

where the indices are taken modulo $n$, i.e., $\mathbf{x}_{i+n} = \mathbf{x}_i$.

## Implementation

### Function `circulant_kronecker`

```@repl
using UniqueKronecker
x = [1, 2];
y = [3, 4];
z = [5, 6];
circulant_kronecker(x, y, z)
```

```@docs
circulant_kronecker
```

### Operator `⊛`

#### Example 1: Two Vectors

```@repl
using UniqueKronecker
x = [1, 2]
y = [3, 4]
x ⊛ y
```

Explanation:

- The circulant Kronecker product for two vectors is:
  ```math
  \mathbf{x} \circledast \mathbf{y} = \mathbf{x} \otimes \mathbf{y} + \mathbf{y} \otimes \mathbf{x}
  ```
- Computation:
  ```julia
  x ⊗ y = [1*3, 1*4, 2*3, 2*4] = [3, 4, 6, 8]
  y ⊗ x = [3*1, 3*2, 4*1, 4*2] = [3, 6, 4, 8]
  Sum: [3+3, 4+6, 6+4, 8+8] = [6, 10, 10, 16]
  ```

Note: The output `[6, 10, 10, 16]` is the sum of the Kronecker products in this case.

#### Example 2: Three Vectors

```@repl
using UniqueKronecker
x = [1, 2]
y = [3, 4]
z = [5, 6]
⊛(x, y, z)
```

Explanation:

- The circulant Kronecker product for three vectors is:
  ```math
  \mathbf{x} \circledast \mathbf{y} \circledast \mathbf{z} = \mathbf{x} \otimes \mathbf{y} \otimes \mathbf{z} + \mathbf{y} \otimes \mathbf{z} \otimes \mathbf{x} + \mathbf{z} \otimes \mathbf{x} \otimes \mathbf{y}
  ```
- Computation involves computing each Kronecker product and summing the results.

#### Example 3: Matrices

```@repl
using UniqueKronecker
A = [1 2; 3 4];
B = [5 6; 7 8];
C = [9 10; 11 12];
result = ⊛(A, B, C)
```

## Generalization

With this implementation, you can now use the ⊛ operator or `circulant_kronecker` function with any number of vectors or matrices. The operation will automatically handle the necessary cyclic permutations and compute the sum of their Kronecker products.

## Columnwise Operation on Snapshot Matrices

Using the function `circulant_kron_snapshot_matrix`, the circulant Kronecker product operation can be applied to construct a snapshot matrix with each column containing the circulant Kronecker product result computed from each column of the input matrices. 

### Example

```@repl
using UniqueKronecker
X1 = [1 2; 3 4]
X2 = [5 6; 7 8]
X3 = [9 10; 11 12]
X = circulant_kron_snapshot_matrix(X1, X2, X3)
```

```@docs
circulant_kron_snapshot_matrix
```

## Derivative of Kronecker products

This circulant Kronecker product is useful because it is a well-known result that the derivative of Kronecker products hold the circulant Kronecker structure. For example, if $\mathbf{x}\in\mathbb{R}^2$ then

```math
\frac{\partial}{\partial \mathbf{x}}(\mathbf{x} \otimes \mathbf{x}) = \mathbf{I}_2 \otimes \mathbf{x} + \mathbf{x} \otimes \mathbf{I}_2
```

and 

```math
\frac{\partial}{\partial \mathbf{x}}(\mathbf{x} \otimes \mathbf{x} \otimes \mathbf{x}) = \mathbf{I}_2 \otimes \mathbf{x} \otimes \mathbf{x} + \mathbf{x} \otimes \mathbf{I}_2 \otimes \mathbf{x} + \mathbf{x} \otimes \mathbf{x} \otimes \mathbf{I}_2
```

which can be reformulated using the circulant Kronecker product as, e.g.,

```math
\frac{\partial}{\partial \mathbf{x}}(\mathbf{x} \otimes \mathbf{x}) = \mathbf{I}_2 \circledast \mathbf{x}
```

## Symmetry and Conversion using Elimination and Duplication Matrices

The circulant Kronecker product inherently captures the symmetric structure present in Kronecker products involving permutations of vectors or matrices. This symmetry can be leveraged to simplify computations, reduce redundancy, and improve efficiency in mathematical operations, particularly when dealing with higher-order tensors or polynomial systems.

### Role of Elimination and Duplication Matrices

When working with circulant Kronecker products, especially in the context of symmetric tensors or polynomial dynamical systems, you may need to convert between the full Kronecker product and its unique representation.

#### From Full Kronecker Product to Unique Representation

1. **Compute the Full Circulant Kronecker Product**: Use the `circulant_kronecker` function to compute the sum over all cyclic permutations.

2. **Apply the Symmetrizer Matrix**: Multiply the result by the symmetrizer matrix $\mathbf{S}_{n,k}$ to make the redundant elements symmetric. Unlike the unique Kronecker product, this must be explicitly applied for circulant Kronecker products.

3. **Apply the Elimination Matrix**: Multiply the result by the elimination matrix $\mathbf{L}_{n,k}$ to extract the unique elements.

   ```math
   \text{Unique Elements} = \mathbf{L}_{n,k} \times \mathbf{S}_{n,k} \times \left( \mathbf{x}_1 \circledast \mathbf{x}_2 \circledast \dots \circledast \mathbf{x}_n \right)
   ```

#### From Unique Representation to Full Kronecker Product

1. **Start with Unique Elements**: Obtain the vector containing the unique monomials.

2. **Apply the Duplication Matrix**: Multiply the unique elements by the duplication matrix $\mathbf{D}_{n,k}$ to reconstruct the full symmetric tensor.

   ```math
   \text{Full Kronecker Product} = \mathbf{D}_{n,k} \times \text{Unique Elements}
   ```

##### Example

Consider vectors $\mathbf{x}, \mathbf{y} \in \mathbb{R}^2$:

1. **Compute the Circulant Kronecker Product**:

   ```math
   \mathbf{x} \circledast \mathbf{y} = \mathbf{x} \otimes \mathbf{y} + \mathbf{y} \otimes \mathbf{x}
   ```

2. **Apply the Symmetrizer and Elimination Matrix**:

   - Use $\mathbf{S}_{2,2}$ and $\mathbf{L}_{2,2}$ to extract unique terms:

     ```math
     \text{Unique Elements} = \mathbf{L}_{2,2}~\mathbf{S}_{2,2} (\mathbf{x} \circledast \mathbf{y})
     ```

3. **Result**:

   - The unique elements correspond to the monomials $x_1 y_1$, $x_1 y_2 + x_2 y_1$, and $x_2 y_2$.
