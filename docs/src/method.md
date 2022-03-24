# Method
Here, we detail CFMM's algorithmic approach to the optimal arbitrage problem
$$
\begin{array}{ll}
\text{maximize}     & U(\Psi) \\
\text{subject to}   & \Psi = \sum_{i=1}^m A_i(\Lambda_i - \Delta_i) \\
& \phi_i(R_i + \gamma_i\Delta_i - \Lambda_i) \geq \phi_i(R_i), \quad i = 1, \dots, m \\
&\Delta_i \geq 0, \quad \Lambda_i \geq 0, \quad i = 1, \dots, m.
\end{array}
$$

## Partial Dualization
It's easy to see that only the first constraint couples this problem; all other 
constraints are specific to a particular CFMM. Thus, if we get rid of this
constraint, the CFMM subproblems can be solved separately. 
Motivated by this insight, we form the Lagrangian

$$
\mathcal{L}(\Psi, \Delta, \Lambda, \nu) = U(\Psi) +
\nu^T(\Psi - \sum_{i=1}^m A_i(\Lambda_i - \Delta_i)) + 
\sum_{i=1}^m \mathbf{I}_i(\Delta_i, \Lambda_i),
$$
where
$$
\mathbf{I}_i(\Delta_i, \Lambda_i) =
\begin{cases}
0, \qquad \Delta_i, \Lambda_i \geq 0, \text{ and } \phi_i(R_i + \gamma_i\Delta_i - \Lambda_i) \geq \phi(R_i) \\
\infty, \qquad \text{o.w.}
\end{cases}
$$
The corresponding dual function is
$$
\begin{align*}
g(\nu) &= \inf_{\Psi, \Delta_i, \Lambda_i} \mathcal{L}(\Psi, \Delta, \Lambda, \nu) \\
&= \inf_\Psi \left(U(\Psi) +
\nu^T\Psi\right) + \sum_{i=1}^m \inf_{\Delta_i, \Lambda_i}\left(
    \mathbf{I}_i(\Delta_i, \Lambda_i) - (A_i^T\nu)^T(\Lambda_i - \Delta_i) \right) \\
    &= -(U^*(-\nu)) + \sum_{i=1}^m -\mathbf{arb}_i(\nu),
\end{align*}
$$
where $\mathbf{arb}_i(\nu)$ is the optimal arbitrage profit from CFMM $i$ with
external market prices $\nu$. These can be evaluated in parallel by solving the convex problem
$$
\begin{array}{ll}
\text{maximize}     & (A_i^T\nu)^T(\Lambda_i - \Delta_i)\\
\text{subject to}   & \phi_i(R_i + \gamma_i\Delta_i - \Lambda_i) \geq \phi_i(R_i) \\
&\Delta_i \geq 0, \quad \Lambda_i \geq 0.
\end{array}
$$

This problem can be efficiently solved with standard techniques (e.g., a primal-dual iterior point method). 
For some common CFMMS (e.g., Uniswap v2), this problem has a closed form solution [CITE]. 

## Minimizing the Dual Problem
As a result of convex duality theory, minimizing $g(\nu)$ is equivalent to solving the original problem in that the optimal values.
We take this dual approach to take advantage of the problem's separability.

To minimize $g(\nu)$, we use LBFGS-B [CITE], which requires evaluation of $g(\nu)$
and $\nabla g(\nu)$. Note that the gradient the Fenchel conjugate $\nabla f^*(y)$ is the 
$x$ at which the supremum $\sup_x (y^Tx - f(x))$ is attained (this could be a set of points--called the subdifferential--in which case we simply choose any $x$ in this set).

The gradient $\nabla_\nu \mathbf{arb}_i(\nu) = A_i  (\Lambda_i^* - \Delta_i^*)$, where $\Lambda^*$ and $\Delta_i^*$ are the optimal values associated with $\mathbf{arb}_i(\nu)$. Thus,
$$
\begin{align}
\nabla g(\nu) &= \Psi^* + \sum_{i=1}^m -A_i  (\Lambda_i^* - \Delta_i^*).

\end{align}
$$
