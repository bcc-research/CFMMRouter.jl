using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using LinearAlgebra, Random, StatsBase, BenchmarkTools
using Plots, LaTeXStrings
using CFMMRouter

## Setup experiment
ns_pools = round.(Int, 10 .^ range(2, 4, 10))
factors = [1, 2, 4]
times = zeros(length(ns_pools), length(factors))
σs = zeros(length(ns_pools), length(factors))

function run_trial(n_pools, factor; rseed=1234)
    Random.seed!(rseed)
    n_tokens = round(Int, factor * sqrt(n_pools))
    v0 = ones(n_tokens)

    pools = Vector{CFMM{Float64}}(undef, n_pools)
    coins = collect(1:n_tokens)

    for i in 1:length(pools)
        Ai = sample(coins, 2, replace=false)
        pools[i] = ProductTwoCoin(
            1000*rand(2),
            rand((0.997, 1.0)),
            Ai
        )
    end

    router = Router(
        LinearNonnegative(rand(n_tokens)),
        pools,
        n_tokens
    )

    trial = @benchmark route!($router; v=$v0)
    return trial
end

## Run experiment
for (i, n_pools) in enumerate(ns_pools)
    for (j, factor) in enumerate(factors)
        trial = run_trial(n_pools, factor)
        
        times[i, j] = median(trial).time
        σs[i, j] = std(trial).time
    end
    @info "Done with $n_pools"
end

## Plot results
plt = plot(
    ns_pools,
    times ./ 1e9,
    lw=3,
    ribbon=σs ./ 1e9,
    yaxis=:log,
    xaxis=:log,
    fillalpha=0.5,
    title="Routing Solve Time",
    ylabel="Time (seconds)",
    xlabel="Number of Swap Pools (m)",
    legend=:bottomright,
    label=[L"$\sqrt{m}$ tokens" L"$2\sqrt{m}$ tokens" L"$4\sqrt{m}$ tokens"],
    dpi=300,
    xlims=(100, 10_000),
    minorgrid=true,
    margin=3Plots.PlotMeasures.mm,
)
savefig(plt, joinpath(@__DIR__, "router_scaling.png"))