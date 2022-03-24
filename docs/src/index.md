```@meta
CurrentModule = CFMMRouter
```

# CFMM Router
Convex optimization for fun and profit.

```@contents
```

## Overview

`CFMMRouter.jl` provides a high performance implementation of the optimal multi-exchange routing problem 
from [Optimal Routing for Constant Function Market Makers](https://web.stanford.edu/~guillean/papers/cfmm-routing.pdf) [1]:

```math
\begin{array}{ll}
\text{maximize}     & U(\Psi) \\
\text{subject to}   & \Psi = \sum_{i=1}^m A_i(\Lambda_i - \Delta_i) \\
& \phi_i(R_i + \gamma_i\Delta_i - \Lambda_i) = \phi_i(R_i), \quad i = 1, \dots, m \\
&\Delta_i \geq 0, \quad \Lambda_i \geq 0, \quad i = 1, \dots, m.
\end{array}
```
A routing problem is specified by the list of CFMMs 
(each with an associated subset of tokens, trading function $\phi_i$, reserves $R_i$ and fee $\gamma_i$), 
and the concave utility function $U(\Psi)$.
We solve for the tendered baskets $\Delta_i \in \mathbf{R}^{n_i}_+$, the received baskets $\Lambda_i \in \mathbf{R}^{n_i}_+$, and the net trade $\Psi \in \mathbf{R}^n$.

Once this solution is obtained, a list of trades can be constructed with
a small amount of additional effort. 
We can generate this list by solving a network flow problem with lossy edges,
but we'll leave the details for a future post.


## Optimization Approach
Our approach is a common one in large scale optimization: we decompose the problem such that the optimal trade on each individual CFMM can be solved independently. It is easy to see that the only constraint coupling the CFMM trades together is

```math
\Psi = \sum_{i=1}^m A_i(\Lambda_i - \Delta_i).
```
Therefore, if we can eliminate this constraint, the subproblems can be solved independently. Intuitively, our approach is to relax this constraint to be a penalty in the objective, where there is some cost of violation. If we fix these costs at a given round, the CFMM subproblems can be solved independently. We then use these subproblem solutions to update the cost of violation and iterate this process. The algorithmic details are explained in [Method](@ref).

### References
[1] G Angeris, T Chitra, A Evans, and S Boyd. [Optimal Routing for Constant Function Market Makers](https://web.stanford.edu/~guillean/papers/cfmm-routing.pdf)