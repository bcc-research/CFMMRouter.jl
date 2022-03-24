# Method
Here we describe our algorithmic approach to solve the routing problem
```math
\begin{array}{ll}
\text{maximize}     & U(\Psi) \\
\text{subject to}   & \Psi = \sum_{i=1}^m A_i(\Lambda_i - \Delta_i) \\
& \phi_i(R_i + \gamma_i\Delta_i - \Lambda_i) = \phi_i(R_i), \quad i = 1, \dots, m \\
&\Delta_i \geq 0, \quad \Lambda_i \geq 0, \quad i = 1, \dots, m.
\end{array}
```

## Overview
Our approach is a common one in large scale optimization: we decompose the problem such that the optimal trade on each individual CFMM can be solved independently. It is easy to see that the only constraint coupling the CFMM trades together is
```math
\Psi = \sum_{i=1}^m A_i(\Lambda_i - \Delta_i).
```
Therefore, if we can eliminate this constraint, the subproblems can be solved independently. Intuitively, our approach is to relax this constraint to be a penalty in the objective, where there is some cost of violation. If we fix these costs at a given round, the CFMM subproblems can be solved independently. We then use these subproblem solutions to update the cost of violation and iterate this process. The algorithmic details are explained below.

## Partial Dualization of the Constraints
Motivated by the insight above, we form the (partial) Lagrangian

```math
\mathcal{L}(\Psi, \Delta, \Lambda, \nu) = U(\Psi) +
\nu^T(\Psi - \sum_{i=1}^m A_i(\Lambda_i - \Delta_i)) +
\sum_{i=1}^m \mathbf{I}_i(\Delta_i, \Lambda_i),
```
where
```math
\mathbf{I}_i(\Delta_i, \Lambda_i) =
\begin{cases}
0, \qquad \Delta_i, \Lambda_i \geq 0, \text{ and } \phi_i(R_i + \gamma_i\Delta_i - \Lambda_i) \geq \phi(R_i) \\
\infty, \qquad \text{o.w.}
\end{cases}
```
Notice that we only introduce a dual variable for the coupling constraint. The corresponding dual function, found by minimization over the primal variables, is

```math
\begin{aligned}
g(\nu) &= \inf_{\Psi, \Delta_i, \Lambda_i} \mathcal{L}(\Psi, \Delta, \Lambda, \nu) \\
&= \inf_\Psi \left(U(\Psi) +
\nu^T\Psi\right) + \sum_{i=1}^m \inf_{\Delta_i, \Lambda_i}\left(
\mathbf{I}_i(\Delta_i, \Lambda_i) - (A_i^T\nu)^T(\Lambda_i - \Delta_i) \right) \\
&= -(U^*(-\nu)) + \sum_{i=1}^m -\mathbf{arb}_i(\nu),
\end{aligned}
```

where $\mathbf{arb}_i(\nu)$ is the optimal arbitrage profit from CFMM $i$ with external "market prices" $\nu$ and $U^*$ is the Fenchel conjugate of $U$. These can be evaluated in parallel by solving the convex problem

```math
\begin{array}{ll}
\text{maximize} & (A_i^T\nu)^T(\Lambda_i - \Delta_i)\\
\text{subject to} & \phi_i(R_i + \gamma_i\Delta_i - \Lambda_i) \geq \phi_i(R_i) \\
&\Delta_i \geq 0, \quad \Lambda_i \geq 0.
\end{array}
```

This arbitrage problem can be efficiently solved with standard techniques (e.g., a primal-dual interior point method). Furthermore, for many common CFMMS (e.g., Uniswap v2), this problem has a closed form solution (for example, see [2, Appendix A]), which we exploit in our solver.

## Minimizing the Dual Problem
As a result of strong duality, minimizing $g(\nu)$ is equivalent to solving the original problem in that the optimal values are equal. After optimizing $g(\nu)$, we can reconstruct the trade using the $\Delta_i$'s and $\Lambda_i$'s found in the sub problems.

To minimize $g(\nu)$, we use LBFGS-B [3], which requires evaluation of $g(\nu)$ and $\nabla g(\nu)$. Note that the gradient of the Fenchel conjugate $\nabla f^*(y)$ is the $x$ at which the supremum $\sup_x (y^Tx - f(x))$ is attained (this may not be a unique point but instead be a set of points--called the subdifferential--in which case we simply choose any $x$ in this set).

The gradient $\nabla_\nu \mathbf{arb}_i(\nu) = A_i (\Lambda_i^* - \Delta_i^*)$, where $\Lambda^*$ and $\Delta_i^*$ are the optimal values associated with $\mathbf{arb}_i(\nu)$. Thus,

```math
\begin{aligned}
\nabla g(\nu) &= \Psi^* + \sum_{i=1}^m -A_i (\Lambda_i^* - \Delta_i^*).
\end{aligned}
```

By evaluating the function $g(\nu)$, we get the gradient essentially for free. 

In the near future, we will support user-created CFMMs as well, which are specified by the invariant function $\phi$, it's gradient $\nabla \phi$, and its Hessian $\nabla^2\phi$. 
The gradient and Hessian may be specified exactly or by using automatic differentiation tools such as `ForwardDiff.jl` [4].

## References
1. G Angeris, T Chitra, A Evans, S Boyd (2021). [Optimal Routing for Constant Function Market Makers.](https://web.stanford.edu/~guillean/papers/cfmm-routing.pdf)
2. G Angeris, Kao, T H, R Chiang, C Noyes, T Chitra (2019). [An analysis of Uniswap markets.](https://web.stanford.edu/~guillean/papers/uniswap_analysis.pdf)
3. C Zhu, Byrd, H R, P Lu, J Nocedal (1997). [L-BFGS-B](http://users.iems.northwestern.edu/~nocedal/lbfgsb.html)
4. J Revels, M Lubin, T Papamarkou (2016). [Forward-mode automatic differentiation in Julia.](https://arxiv.org/abs/1607.07892)