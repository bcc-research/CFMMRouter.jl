# Credit: Chris Rackauckas
macro def(name, definition)
    return quote
        macro $(esc(name))()
            esc($(Expr(:quote, definition)))
        end
    end
end

# Adapted from
# https://github.com/ZIB-IOL/FrankWolfe.jl/blob/master/src/fw_algorithms.jl
function print_header(data)
    @printf(
        "\n─────────────────────────────────────────────────────────────────────────────────────────────────\n"
    )
    @printf(
        "%13s %14s %14s %14s %14s %14s\n",
        data[1],
        data[2],
        data[3],
        data[4],
        data[5],
        data[6]
    )
    @printf(
        "─────────────────────────────────────────────────────────────────────────────────────────────────\n"
    )
end

function print_footer()
    @printf(
        "─────────────────────────────────────────────────────────────────────────────────────────────────\n\n"
    )
end

function print_iter_func(data)
    @printf(
        "%13s %14e %14e %14e %14.3e %13.3f\n",
        data[1],
        Float64(data[2]),
        Float64(data[3]),
        Float64(data[4]),
        data[5],
        data[6]
    )
end