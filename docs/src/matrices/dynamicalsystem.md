# Polynomial Dynamical System

The main use case of the unique Kronecker product is in a polynomial structured dynamical systems defined by

```math
    \dot{\mathbf{x}}(t) = \mathbf{A}\mathbf{x}(t) + \mathbf{A}_2\mathbf{x}^{[2]}(t) + \mathbf{A}_3\mathbf{x}^{[3]}(t) + \cdots 
```

where $\mathbf{x}$ is the system variable and $\mathbf{A}_2\in\mathbb{R}^{n\times n^2}$, $\mathbf{A}_3\in\mathbb{R}^{n^3}$ are the system matrices/operators. This can be represented using the unique Kronecker product as 

```math
    \dot{\mathbf{x}}(t) = \mathbf{A}\mathbf{x}(t) + \mathbf{A}_{2u}\mathbf{x}^{\langle 2 \rangle}(t) + \mathbf{A}_{3u}\mathbf{x}^{\langle 3 \rangle}(t) + \cdots  
```

where $\mathbf{A}_{2u}\in\mathbb{R}^{n(n+1)/2}$ and $\mathbf{A}_{3u}\in\mathbb{R}^{n(n+1)(n+2)/6}$. This package is also dedicated to converting between the unique operators $\mathbf{A}_{2u},~\mathbf{A}_{3u}$ and non-unique operators $\mathbf{A}_2,~\mathbf{A}_3$. 


# Creating the Polynomial Operators
This packages offers a built-in function that allows the construction of redundant and non-redundant polynomial operators generalized for the dimension of the system variable and order of the Kronecker product. This function is called `makePolyOp` which easily creates the matrices by accepting the indices and values as arguments.

```@docs
makePolyOp
```

# Example

Here we show some famous systems which possess the polynomial structure:

## Van der Pol Oscillator

The Van der Pol oscillator is a nonlinear system that exhibits self-sustained oscillations with non-constant amplitude and frequency. It is governed by the second-order differential equation:

```math
    \ddot{x} - \mu (1 - x^2) \dot{x} + x = 0
```

where $x$ is the state variable and $\mu$ is a scalar parameter representing the nonlinearity and damping of the system.

To express this as a first-order system suitable for the polynomial dynamical system framework, we define the state vector $\mathbf{x} = [x_1, x_2]^T$, where $x_1 = x$ and $x_2 = \dot{x}$. The system then becomes:

```math
\begin{aligned}
    \dot{x}_1 &= x_2 \\
    \dot{x}_2 &= \mu (1 - x_1^2) x_2 - x_1
\end{aligned}
```

This can be rewritten in matrix form as:

```math
    \dot{\mathbf{x}} = \mathbf{A} \mathbf{x} + \mathbf{A}_{3u} \mathbf{x}^{\langle 3 \rangle}
```

Here, the matrices are defined as:

- Linear term:

  ```math
      \mathbf{A} = \begin{bmatrix}
          0 & 1 \\
          -1 & \mu
      \end{bmatrix}
  ```

- Cubic term:

  ```math
      \mathbf{A}_{3u} = \begin{bmatrix}
          0 & 0 & 0 & 0 \\
          0 & -\mu & 0 & 0
      \end{bmatrix}
  ```

The cubic term arises from the $- \mu x_1^2 x_2$ component, which can be represented using the unique Kronecker product:

```math
    \mathbf{x}^{\langle 3 \rangle} = [x_1^3,~x_1^2 x_2, ~x_1 x_2^2, ~x_2^3]^\top
```

Thus, the Van der Pol oscillator fits into the polynomial dynamical system framework by capturing its nonlinear dynamics through the unique Kronecker product, facilitating analysis and simulation.

!!! tip 
    `UniqueKronecker.jl` provides a convenient function to construct the polynomial operators called `makePolyOp`. For example, you can construct the (redundant) cubic operator using the following command:

    ```@repl
    using UniqueKronecker: makePolyOp
    μ = 1.0
    A3 = makePolyOp(2, [(1,1,2,2)], [μ])
    ```

## Lorenz System

The Lorenz system is a set of three coupled, nonlinear differential equations originally developed to model atmospheric convection. It is famous for its chaotic solutions for certain parameter values. The system is given by:

```math
\begin{aligned}
    \dot{x}_1 &= \sigma (x_2 - x_1) \\
    \dot{x}_2 &= x_1 (\rho - x_3) - x_2 \\
    \dot{x}_3 &= x_1 x_2 - \beta x_3
\end{aligned}
```

where $x_1$, $x_2$, and $x_3$ are the state variables, and $\sigma$, $\rho$, and $\beta$ are parameters.

To express the Lorenz system within the polynomial dynamical system framework, we identify the nonlinear terms and represent them using the unique Kronecker product:

- Quadratic terms: $x_1 x_2$, $x_1 x_3$

We define the state vector $\mathbf{x} = [x_1, x_2, x_3]^\top$ and compute the unique Kronecker powers:

```math
    \mathbf{x}^{\langle 2 \rangle} = [x_1^2, ~x_1 x_2, ~x_1 x_3, ~x_2^2, ~x_2 x_3, ~x_3^2]^\top
```

The system can then be represented as:

```math
    \dot{\mathbf{x}} = \mathbf{A} \mathbf{x} + \mathbf{A}_{2u} \mathbf{x}^{\langle 2 \rangle}
```

Where:

- Linear term:

  ```math
      \mathbf{A} = \begin{bmatrix}
          -\sigma & \sigma & 0 \\
          \rho & -1 & 0 \\
          0 & 0 & -\beta
      \end{bmatrix}
  ```

- Quadratic term:

  ```math
      \mathbf{A}_{2u} = \begin{bmatrix}
          0 & 0 & 0 & 0 & 0 & 0 \\
          0 & 0 & -1 & 0 & 0 & 0 \\
          0 & 1 & 0 & 0 & 0 & 0
      \end{bmatrix}
  ```

In this form, the Lorenz system's nonlinear dynamics are captured using the unique Kronecker product, making it amenable to analysis using polynomial system techniques.

!!! tip
    The quadratic operator can be constructed using the `makePolyOp` function as follows:

    ```@repl
    using UniqueKronecker: makePolyOp
    A2 = makePolyOp(3, [(1,3,2), (1,2,3)], [-1, 1])
    ```

## Viscous Burgers' Equation

The viscous Burgers' equation is a fundamental partial differential equation from fluid mechanics and nonlinear acoustics. It models various physical processes like gas dynamics and traffic flow. The equation is:

```math
    \frac{\partial u}{\partial t} + u \frac{\partial u}{\partial x} = \nu \frac{\partial^2 u}{\partial x^2}
```

where $u(x, t)$ is the velocity field, and $\nu$ is the viscosity coefficient.

To convert this PDE into a finite-dimensional polynomial dynamical system, we discretize the spatial domain using methods like finite differences. Let $\mathbf{u}$ be the vector of discretized values of $u(x, t)$. Assuming periodic boundary condition, the discretized equation becomes:

```math
    \dot{\mathbf{u}} = \mathbf{A} \mathbf{u} + \mathbf{A}_2 (\mathbf{u} \otimes \mathbf{u})
```

Here:

- the discretized Laplacian operator (diffusion term) is repsented by $\mathbf{A}$.
- and the discretized convection operator is $\mathbf{A}_2$ .

By rearranging, we can express the nonlinear term using the unique Kronecker product:

```math
    \dot{\mathbf{u}} = \mathbf{A} \mathbf{u} + \mathbf{A}_{2u} \mathbf{u}^{\langle 2 \rangle}
```

where:

```math
    \mathbf{u}^{\langle 2 \rangle} = [u_1^2, ~u_1 u_2, ~\dots, ~u_n^2]^\top
```

This formulation allows us to apply polynomial system analysis tools to the viscous Burgers' equation, enabling efficient simulation and control design for systems modeled by this equation.

!!! Note
    There are additional functions `makePolyOp_faster` and `makePolyOp_parallel` which are experimental functions designed to speed up the construction of polynomial operators. I am planning to improve these functions and do some benchmarks.


```@docs
makePolyOp_faster
makePolyOp_parallel
```