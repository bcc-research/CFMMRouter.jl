using CFMMRouter
using Documenter
using Literate

# For Plots.jl
# https://discourse.julialang.org/t/plotting-errors-when-building-documentation-using-plots-jl-and-documenter-jl/67849
# ENV["GKSwstype"]="100"

EXCLUDED_EXAMPLES = []

# utility function from https://github.com/JuliaOpt/Convex.jl/blob/master/docs/make.jl
fix_math_md(content) = replace(content, r"\$\$(.*?)\$\$"s => s"```math\1```")

# utility functions from https://github.com/oxfordcontrol/COSMO.jl/blob/master/docs/make.jl
fix_suffix(filename) = replace(filename, ".jl" => ".md")
function postprocess(content)
      """
      The source files for all examples can be found in [/examples](https://github.com/bcc-research/CFMMRouter.jl/tree/main/examples).
      """ * content
end

examples_path = joinpath(@__DIR__, "../examples/")
examples = filter(x -> endswith(x, ".jl") && !in(x, EXCLUDED_EXAMPLES), readdir(examples_path))
build_path =  joinpath(@__DIR__, "src", "examples/")

for example in examples
      Literate.markdown(
        examples_path * example, build_path;
        preprocess = fix_math_md,
        postprocess = postprocess,
        flavor = Literate.DocumenterFlavor(),
        credit = true
    )
end

examples_nav = fix_suffix.(joinpath.("examples", examples))

makedocs(;
    modules=[CFMMRouter],
    authors="Guillermo Angeris, Theo Diamandis",
    repo="https://github.com/bcc-research/CFMMRouter.jl/blob/{commit}{path}#L{line}",
    sitename="CFMMRouter.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://bcc-research.github.io/CFMMRouter.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Solution method" => "method.md",
        "Specifying objectives" => "objective.md",
        "Examples" => examples_nav,
        "API reference" => "api.md"
    ],
)

deploydocs(;
    repo="github.com/bcc-research/CFMMRouter.jl",
    devbranch = "main"
)