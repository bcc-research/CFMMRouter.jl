```@meta
CurrentModule = CFMMRouter
```

# CFMM Router
Convex optimization for fun and profit. (Now in Julia!)

### Documentation Contents:
```@contents
Pages = ["index.md", "method.md", "objective.md"]
Depth = 1
```
##### Examples:
```@contents
Pages = ["examples/arbitrage.md", "examples/liquidate.md"]
Depth = 1
```

## Overview
`CFMMRouting.jl` provides a fast, efficient solver for solving routing problems across decentralized exchanges that are implemented as [*constant function market makers*](https://web.stanford.edu/~guillean/papers/cfmm-chapter.pdf) (CFMMs, for short) such as Uniswap, Balancer, or Curve.

## Getting started
To get a quick start, we recommend jumping over to one of the examples, such as: [Liquidating a basket of tokens](@ref).

## Why do we need routing?
In general, decentralized exchanges, such as Uniswap, Balancer, Curve, among many others, are organized as CFMMs—a particular type of automated market maker, whose behavior is mathematically defined by a *trading function*. A trader can then use a CFMM to trade a basket of assets (sometimes known as tokens) for another desired basket of assets, so long as the desired trade satisfies some conditions.

While there may be many assets in a network, in practice, individual CFMMs trade only a small number of these assets. This fact leads to several interesting problems when trading.

In one simple scenario, a trader might want to trade some amount $x$ of token A for token B. It is then reasonable to ask: what is the largest amount of token B that can be received given this amount $x$ of token A? If the trader is only allowed to trade with a single market of their choosing, the solution to this problem is easy—the trader simply checks how much of token B, given quantity $x$ of token A, they expect to receive from each market that trades tokens A and B,
and picks the market that gives the maximum amount of token B. When they are allowed to trade against any number of markets, as is usually the case in decentralized finance, the problem becomes substantially more complicated. 
In this case, an optimal solution could involve not only splitting an order across numerous individual markets trading tokens A and B, but also including sequential trades across any combination of markets which include other tokens.
Finding the best way to execute this order is known as the [*optimal routing problem*](https://angeris.github.io/papers/cfmm-routing.pdf).

The main idea behind this package is to note that the problem of optimal routing can often be phrased as a [convex optimization problem](https://www.stanford.edu/~boyd/cvxbook/). Problems of this form can generally be quickly and robustly solved in practice, even for relatively large problem instances. This package also exploits some additional structure present in the optimal routing problem to further speed up problem solving. More details about this approach can be found in the sections below.

While the solution method is somewhat technical, using this package requires very little optimization knowledge. We highly recommend starting with the examples provided!

## Technical overview
This package provides a high performance implementation of the optimal multi-exchange routing problem 
from [Optimal routing for constant function market makers](https://web.stanford.edu/~guillean/papers/cfmm-routing.pdf)[^1]:

```math
\begin{array}{ll}
\text{maximize}     & U(\Psi) \\
\text{subject to}   & \Psi = \sum_{i=1}^m A_i(\Lambda_i - \Delta_i) \\
& \phi_i(R_i + \gamma_i\Delta_i - \Lambda_i) = \phi_i(R_i), \quad i = 1, \dots, m \\
&\Delta_i \geq 0, \quad \Lambda_i \geq 0, \quad i = 1, \dots, m.
\end{array}
```
A routing problem is specified by the list of CFMMs that the trader can interact with and a concave, increasing utility function $U: \mathbf{R}^n \to \mathbf{R}$, depending only on the *net output* of all trades (we will define what this means momentarily).

Each CFMM is associated with a subset of $n_i$ tokens which we will write as some subset $S_i \subseteq \{1, \dots, n\}$. Every CFMM $i$ has some trading function $\phi_i: \mathbf{R}^{n_i} \to \mathbf{R} \cup \{\infty\}$, reserves $R_i \in \mathbf{R}^{n_i}$, fee $0 < \gamma_i \le 1$, and matrix $A_i \in \mathbf{R}^{n \times n_i}$, mapping the local token basket that the CFMM trades to the network's token basket. We solve for the tendered baskets $\Delta_i \in \mathbf{R}^{n_i}_+$, the received baskets $\Lambda_i \in \mathbf{R}^{n_i}_+$ for each CFMM $i$, which results in the network trade vector $\Psi \in \mathbf{R}^n$ that maximizes the function $U$. For more information on this problem, we recommend reading the original paper[^1].

#### Uses
The routing problem includes a large number of tasks that users might be interested in, including, for example
- Liquidating a basket of tokens
- Trading token A for token B
- Finding arbitrage opportunities
- And more...

#### Solutions and net trades
A *solution* for this problem is provided as an (unordered) set of trades
$(\Delta_i, \Lambda_i)$ to be performed with each CFMM $i$.

If *all* trades provided by the solution were to be executed, some tokens might have negative balances (*i.e.*, the user needs to tender these tokens for the transaction to succeed) or positive balances (*i.e.*, the user will instead *receive* this amount of tokens, assuming the trades all succeed). The resulting amounts of each token to be tendered or received, after all trades are netted out, is called the *network trade vector*, which we have denoted $\Psi$ above.

#### Executing solutions
The solution this package provides are a set of trades that, when netted out, maximizes the utility function $U$.  With a small amount of additional effort, it is possible to convert this unordered set of trades into a list of trades that can be included in a transaction. We will cover some approaches for how to convert a set of trades into a reasonable transaction in a future post.


## Optimization approach
To solve the optimal routing problem presented here, this package uses a common technique in large scale optimization: we decompose the problem such that the optimal trade on each individual CFMM can be solved independently. It is easy to see that the only constraint coupling the CFMM trades together is

```math
\Psi = \sum_{i=1}^m A_i(\Lambda_i - \Delta_i).
```
Therefore, if we can eliminate this constraint, the subproblems can be solved in a fully parallelizable way. Intuitively, our approach is to relax this constraint to be a penalty in the objective, where there is some cost of violation. If we fix these costs at a given round, the CFMM subproblems can be solved independently. We then use these subproblem solutions to update the cost of violation and iterate this process. The algorithmic details are explained in the [Solution method](@ref) section.

## References
[^1]: G. Angeris, T. Chitra, A. Evans, and S. Boyd. [Optimal routing for constant function market makers](https://angeris.github.io/papers/cfmm-routing.pdf).