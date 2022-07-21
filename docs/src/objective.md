# Objectives
We include a few simple utility functions with this package. 
Arbitrary utility functions (which can include constraints) are specified in
the 
[objectives.jl](https://github.com/bcc-research/CFMMRouter.jl/blob/main/src/objectives.jl) file.
We outline the included utility functions and how to specify new ones below.

## Included utility functions

### LinearNonnegative (arbitrage)
The [`LinearNonnegative`](@ref) objective is
```math
    U(\Psi) = c^T\Psi - \mathbf{I}(\Psi \geq 0),
```
where $c$ is a positive price vector. $\mathbb{I}(\Psi \ge 0)$ is an indicator function that is $0$ if $\Psi \ge 0$ and $+\infty$ otherwise. This objective requires no net input and a linear utility for each token, defined by the price vector. A nonzero solution to this problem finds an arbitrage in the network—a set of trades that yields strictly positive value, but requires no net input.


### BasketLiquidation (liquidations and swaps)
The [`BasketLiquidation`](@ref) objective is 
```math
    U(\Psi) = \Psi_i - \mathbf{I}(\Psi_{-i} + Δ^\mathrm{in}_{-i} = 0, ~ \Psi_i \geq 0),
```
where $Δ^\mathrm{in}$ is a basket of tokens to be liquidated into token $i$. Here, $\Psi_{-i}$ is the vector $\Psi$ with entry $i$ removed, and $\mathbf{I}(\Psi_{-i} + Δ^\mathrm{in}_{-i} = 0, ~\Psi_i \geq 0)$ is the indicator function which is zero if all conditions are met and is $\infty$, otherwise. In the special case where $Δ^\mathrm{in}_k$ is zero at all indices except for some index $j \ne i$, this objective defines a swap from token $j$ to token $i$ where we attempt to maximize the amount of token $i$ received.
We implement a [`Swap`](@ref) objective as shorthand for this special case.


## Specifying new utility functions
New utility functions can be easily specified following the examples in [objectives.jl](https://github.com/bcc-research/CFMMRouter.jl/blob/main/src/objectives.jl).
We only need to specify the conjugate $f(\nu)$ and its gradient:
```math
    f(\nu) = \sup_\Psi \left(U(\Psi) - \nu^T \Psi \right)
```
and
```math
\nabla f(\nu) = -\Psi^\star,
```
where $\Psi^\star$ is the optimal value of the supremum in $f(\nu)$. In addition, we specify lower bound `lower_limit` and upper bound `upper_limit` for the objective.

Sometimes $f(\nu)$ does not have a closed form solution, so we need to solve an additional optimization problem to evaluate the objective.  For example, consider Markowitz portfolio rebalancing:
```math
U(\Psi) = \mu^T(\Psi + \Delta^\mathrm{in}) - \frac{\gamma}{2}
(\Psi + \Delta^\mathrm{in})^T\Sigma(\Psi + \Delta^\mathrm{in}) 
- \mathbb{I}(\Psi + \Delta^\mathrm{in} \geq 0),
```
where $\mu$ and $\Sigma$ are the (estimated) mean and covariance returns for each token.