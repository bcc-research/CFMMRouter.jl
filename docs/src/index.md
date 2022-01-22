```@meta
CurrentModule = CFMMRouter
```

# CFMM Router
Convex optimization for fun and profit.


## Overview

`CFMMRouter.jl` provides a high performance implementation of the optimal multi-exchange arbitrage from [Optimal Routing for Constant Function Market Makers](https://web.stanford.edu/~guillean/papers/cfmm-routing.pdf) [1]:

$$
\begin{array}{ll}
\text{maximize}     & U(\Psi) \\
\text{subject to}   & \Psi = \sum_{i=1}^m A_i(\Lambda_i - \Delta_i) \\
& \phi_i(R_i + \gamma_i\Delta_i - \Lambda_i) = \phi_i(R_i), \quad i = 1, \dots, m \\
&\Delta_i \geq 0, \quad \Lambda_i \geq 0, \quad i = 1, \dots, m.
\end{array}
$$
An arbitrage problem is specified by the CFMMs $\phi_i$ (with associates reserves and fees),
the coins in CFMM $i$ (captured by $A_i$), and the concave utility function $U(\Psi)$.
We solve for the tendered baskets $\Delta_i \in \mathbf{R}^{n_i}_+$, the received baskets $\Lambda_i \in \mathbf{R}^{n_i}_+$, and the net trade $\Psi \in \mathbf{R}^n$.









## References
[1] Angeris, Chitra, Evans, and Boyd. [Optimal Routing for Constant Function Market Makers](https://web.stanford.edu/~guillean/papers/cfmm-routing.pdf)