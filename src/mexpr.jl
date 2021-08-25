#   This file is part of Maxima.jl. It is licensed under the MIT license
#   Copyright (c) 2016 Nathan Smith
#   Copyright (c) 2021 Jingyue Wang
#
#   TODO: type-check. for example, if we define `b : matrix([1,2],[3,4])`, then
#   `cauchy_matrix([a,b,c])` gives very wierd result that can not be rendered
#   correctly as LaTeX on DomTerm.
export
    MExpr,
    @mx_str,
    @jx_str,
    parse,
    convert,
    get_mstr,
    clean_matrix_env,
    mexpr_to_tex,
    jexpr_to_mstr

import
    Base: parse, convert, *

using MathMLRepr

const infix_ops = [:+, :-, :*, :/, :^, ://, :.*, :∘]

isinfix(args) = args[1] in infix_ops && length(args) > 2

const jsym_to_m = Dict(
    :α           => "alpha",
    :β           => "beta",
    :δ           => "delta",
    :ϵ           => "epsilon",
    :ζ           => "zeta",
    :η           => "eta",
    :θ           => "theta",
    :ι           => "iota",
    :κ           => "kappa",
    :λ           => "lambda",
    :μ           => "mu",
    :ν           => "nu",
    :ξ           => "xi",
    :ρ           => "rho",
    :σ           => "sigma",
    :τ           => "tau",
    :υ           => "upsilon",
    :χ           => "chi",
    :ψ           => "psi",
    :ω           => "omega",
    :ℯ           => "%e",
    :pi          => "%pi",
    :π           => "%pi",
    :γ           => "%gamma",
    :eulergamma  => "%gamma",
    :golden      => "%phi",
    :φ           => "%phi",
    :im          => "%i",
    :Inf         => "inf"
)

jexpr_to_mstr(x::Symbol) = x in keys(jsym_to_m) ? jsym_to_m[x] : string(x)
jexpr_to_mstr(x::Rational) = string(x.num)*"/"*string(x.den)
jexpr_to_mstr(x::Union{AbstractFloat, Integer}) = string(x)
jexpr_to_mstr(x::String) = x
jexpr_to_mstr(x::LineNumberNode) = ""
jexpr_to_mstr(x::QuoteNode) = error("Can not process QuoteNode")
# return a string containing a valid Maxima statement
# We assume expr is NOT a block expression, so basically just one statement.
# The function can handle block expressions inside a function object.
function jexpr_to_mstr(expr::Expr; defeqn::Bool=false)
    o = ""
    if expr.head == :call
        operator = isinfix(expr.args) ? " $(expr.args[1]) " : ", "
        # operators that need special treatments to make a Maxima expression
        if operator == " // "; operator = "/"; end
        if operator == " ∘ "; operator = "."; end
        if operator == " .* "; operator = "*"; end
        if !isinfix(expr.args)
            o *= jexpr_to_mstr(expr.args[1])
        end
        o *= "("* (length(expr.args) > 1 ? jexpr_to_mstr(expr.args[2]) : "")
        for arg in expr.args[3:end] # if s > end args[3:end] is empty
            o *= operator*jexpr_to_mstr(arg)
        end
        o *= ")" #JW: we can replace `()` by `{}` for a component behind `^` in TeX string
    elseif expr.head == :(=)
        #TODO: save variable names and their type info
        if defeqn === true
            o *= jexpr_to_mstr(expr.args[1])*" = "*jexpr_to_mstr(expr.args[2])
        elseif expr.args[1] isa Symbol 
            o *= jexpr_to_mstr(expr.args[1])*": "*jexpr_to_mstr(expr.args[2])
        else # should be a :call that implies a functin definition, e.g., f(x)=x^2
            fun_body = join(jexpr_to_mstrv(expr.args[2]), ",")
            o *= jexpr_to_mstr(expr.args[1])*":= "*fun_body
        end
    elseif expr.head == :block
        # we deal a limited case of :block, namely only one expression inside
        length(expr.args) > 2 && error("The RHS of `=` is too compiated to process")
        # args[1] is assumed to be a LineNumberNode
        !(expr.args[1] isa LineNumberNode) && error("Unknown expression in RHS")
        o *= jexpr_to_mstr(expr.args[2])
    elseif expr.head == :function
        o *= jexpr_to_mstr(expr.args[1])*" := block([]"
        # args[2] of a function always returns a Julia block expression, even for an
        # empty function (`function f(x) end`), it returns an empty block
        args = jexpr_to_mstrv(expr.args[2])
        for a in args; o *= ", "*a; end
        o *= ")"
    elseif expr.head == :row || expr.head == :vect
        o = "["*jexpr_to_mstr(expr.args[1])
        for a in expr.args[2:end]
            o *= ", "*jexpr_to_mstr(a)
        end
        o *= "]"
    elseif expr.head == :vcat
        o = "matrix("
        for a in expr.args[1:end-1]
            o *= jexpr_to_mstr(a)*", "
        end
        o *= jexpr_to_mstr(expr.args[end])*")"
    elseif expr.head == Symbol("'")
        o = "transpose("*jexpr_to_mstr(expr.args[1])*")"
    else
        # can not deal with QuoteNode either
        error("Nested :$(expr.head) block structure not supported by Maxima.jl")
    end

    return o
end


# for a single expression, return an array that only contains one `mstring`
# for a block expression, return an array that contains all expressions in the block
function jexpr_to_mstrv(expr::Expr)
    if expr.head in (:block, :toplevel)
        strvec = String[]
        for a in expr.args
            push!(strvec, jexpr_to_mstr(a))
        end
        # a LineNumeberNode yields an empty string.
        return clean_empty_ele(strvec)
    else
        return String[jexpr_to_mstr(expr)]
    end
end

function jexpr_to_mstr(x::Vector{T}) where T
    o = "["*jexpr_to_mstr(x[1])
    for a in x[2:end]
        o *= ", "*jexpr_to_mstr(a)
    end
    o *= "]"
end

function jexpr_to_mstr(Mat::Matrix{T}) where T
    m, n = size(Mat)
    o = "matrix("
    # add the first row
    o *= jt2ms_addrow(Mat[1, :])
    # add each of the rest of rows
    for i = 2:m
        o *=  ", "*jt2ms_addrow(Mat[i, :])
    end
    o *= ")"
    return o
end

# add a Julia vector as a row to a Maxima table
function jt2ms_addrow(v::Vector{T}) where T
    o = "["*jexpr_to_mstr(v[1])
    for j = 2:length(v)
        o *= ", "*jexpr_to_mstr(v[j])
    end
    o *= "]"
    return o
end

"""

A Maxima expression

## Summary:
    
type MExpr <: Any

## Fields:

strvec :: Vector{String}
"""
struct MExpr
    strvec::Vector{String}
end
MExpr(x::Rational) = MExpr("$(x.num)"*"/"*"$(x.den)")
MExpr(x::Number) = MExpr(string(x))                           
MExpr(x::String) = MExpr(push!(String[], x))
MExpr(x::Symbol) = MExpr(string(x))
"""
    MExpr(expr::Expr)

Convert a Julia expression `expr` to an equivalent Maxima expression

## Examples
```julia
julia> MExpr(:(sin(x*im) + cos(y*φ)))

                           cos(%phi y) + %i sinh(x)

```
"""
MExpr(expr::Expr) = MExpr(jexpr_to_mstrv(expr))


macro mx_str(x) MExpr(x) end
macro jx_str(x) MExpr(Meta.parse(x)) end


*(x::MExpr,  y::MExpr)  = MExpr(append!(String[], x.strvec, y.strvec))
*(x::MExpr,  y::String) = MExpr(append!(String[], x.strvec, [y]))
*(x::String, y::MExpr)  = MExpr(append!(String[x], y.strvec))

get_mstr(m::MExpr, i) = m.strvec[i]

clean_empty_ele(m::Vector{String}) = filter(!isempty, m)
allspace(x::AbstractString) = all(isspace, x)
clean_space_ele(m::Vector{String}) = filter(!allspace, m)
# for each mstring in a MExpr, if it is a command sequence, split it into
# a vector of single-commands. Maxima commands are separated by `$` or `;` in
# a sequence.
function normalize(m::MExpr)
    r = String[]
    for s in m.strvec
        p = split(s, ['$', ';'])
        # NOTE: `split()` returns `p`, a Vector{SubString{String}}, NOT a 
        #Vector{String}. `append!()` makes `r` a Vector{String}
        append!(r, p)
    end
    return MExpr(r |> clean_empty_ele |> clean_space_ele)
end


const m_to_jstr = Dict(
    "alpha"   =>  "α",          
    "beta"    =>  "β",          
    "delta"   =>  "δ",          
    "epsilon" =>  "ϵ",          
    "zeta"    =>  "ζ",          
    "eta"     =>  "η",          
    "theta"   =>  "θ",          
    "iota"    =>  "ι",          
    "kappa"   =>  "κ",          
    "lambda"  =>  "λ",          
    "mu"      =>  "μ",          
    "nu"      =>  "ν",          
    "xi"      =>  "ξ",          
    "rho"     =>  "ρ",          
    "sigma"   =>  "σ",          
    "tau"     =>  "τ",          
    "upsilon" =>  "υ",          
    "chi"     =>  "χ",          
    "psi"     =>  "ψ",          
    "omega"   =>  "ω",          
    "%e"      =>  "ℯ",
    "%pi"     =>  "π",
    "%i"      =>  "im",
    "%gamma"  =>  "eulergamma",
    "%phi"    =>  "φ",
    "inf"     =>  "Inf",
    "minf"    =>  "-Inf"
)
"""
    parse(mexpr::MExpr)

Convert a Maxima expression into an equivalent Julia expression

## Examples
```julia
julia> parse(mx\"sin(%i*x)\")
:(im * sinh(x))

```
"""
#FIXME: I should do real parsing instead of naive text substitutions.
#NOTE: `parse()` lambdify a symbolic function, convert a symbolic number into
#      a concrete type (Int64, Float32,...).
function parse(m::MExpr)
    exprs = Any[]
    for ex in normalize(m).strvec
        for key in keys(m_to_jstr)
            #FIXME: we should parse the MExpr, not simply doing a text substitution.
            ex = replace(ex, key => m_to_jstr[key])
        end
        if occursin( ":=", ex)
            args = split(ex, ":=")
            push!(
                exprs,
                Expr(
                    :function,
                    parse(args[1]),  # e.g. parse("f(x)")
                    args[2] |> String |> MExpr |> parse))
        elseif occursin( ":", ex)
            args = split(ex, ":")
            push!(
                exprs,
                Expr(
                    :(=),
                    parse(args[1]),
                    args[2] |> String |> MExpr |> parse))
        elseif occursin( "block([], ", ex)
            x = replace(ex, "block([],", "") |> strip |> chop  # chop the last `)`
            xs = split(x, ",")
            l = Vector{Any}(length(x))
            for s in xs
                !isempty(s) && push!(l, s |> String |> MExpr |> parse)
            end
            push!(exprs, Expr(:block, l...))
        else
            push!(exprs, Meta.parse(ex))
        end
    end
    return length(exprs) == 1 ? exprs[1] : Expr(:block, exprs...)
end


convert(::Type{MExpr}, m::MExpr) = m
convert(::Type{Array{String,1}}, m::MExpr) = m.strvec
convert(::Type{String}, m::MExpr) = join(m.strvec, "; ")
convert(::Type{T}, m::MExpr) where T<:Number = eval(parse(m))
#IMPORTANT! we strip the leading and traling `$`s generated by Maxima
function convert(::Type{LATEX}, m::MExpr)
    check = "'(" * replace(convert(String, m), r";" => ")\$\n'(") * ")"
    write(ms.input, "$check\$\n print(ascii(4))\$")
    out = (readuntil(ms.output, EOT) |> String
                                     |> str -> rstrip(str, EOT))

    s = ""
    if occursin(synerr, out) || occursin(maxerr, out)
        @warn "Invalid Maxima expression"
        s = "Invalid Maxima Expression: "*out
    else
        write(ms, "tex1('("*replace(convert(String, m), r";" => "))\$\ntex1('(") * "))")
        texstr = read(ms)
        s = replace(texstr, r"\nfalse\n" => "")
        # Warning! quick hacks inside this function
        s = clean_matrix_env(s)
        s = String(strip(s, '$'))
    end
    return LATEX(s)
end
#JW: In stdlib 1.6.1 `display(x::MExpr) dereferences x into Any
convert(::Type{Any}, m::MExpr) = m
mexpr_to_tex(m::MExpr) = convert(LATEX, m)

# All tricks to post-process the TeX code generated by maxima are wrapped
# in this function
function clean_matrix_env(x::String)
    x = strip(x)
    x = chop(x, head=1, tail=1) # x retuend by `tex1()` is encluded by double-quote `"`
    if occursin("pmatrix", x)
        x = replace(x, r"\\\\ifx.*?\\\\fi"=>"") # del all the \ifx
        x = replace(x, r"\\\\cr +$"=>"") # del the last empty \cr
        x = replace(x, r"\\\\cr"=>"\\\\") # \cr -> \\
        x = replace(x, r"\\\\([a-z])"=>s"\\\1") # any double-slash command \\cmd => \cmd
        x = replace(x, r"\\\\([,!])"=>s"\\\1") # any double-slash command \\opr => \opr
        x = "\\begin{pmatrix}"*x*"\\end{pmatrix}"
    else
        x = replace(x, raw"\\\\"=>"\\")
    end
    return x
end
