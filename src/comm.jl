# Communication with Maxima
export 
    MaximaError,
    MaximaSyntaxError,
    mcall,
    mcallj,
    @mx,
    ==

import
    Base: ==



struct MaximaError <: Exception
    errstr::String
end
Base.showerror(io::IO, err::MaximaError) = print(io, err.errstr)

struct MaximaSyntaxError <: Exception
    errstr::String
end
Base.showerror(io::IO, err::MaximaSyntaxError) = print(io, err.errstr)

"""
    mcall(m::MExpr)

Evaluate a Maxima expression.

## Examples
```julia
julia> m\"integrate(sin(x), x)\"

                             integrate(sin(x), x)

julia> mcall(ans)

                                   - cos(x)

```
"""
function mcall(m::MExpr)
    write(ms, convert(String, m))
    output = read(ms)
    if occursin(maxerr, output)
        write(ms.input, "errormsg()\$")
        write(ms.input, "print(ascii(4))\$")
        message = read(ms)
        throw(MaximaError(message))
    elseif occursin(synerr, output)
        throw(MaximaSyntaxError(output))
    else
        r = String[]
        for k in split(output, ETX)
            push!(r, replace(k, r"[ \n]"=>""))
        end
        return normalize(MExpr(r))
    end
end
"""
  mcall(expr::T) where T

Evaluate a Julia expression or string using the Maxima interpretor.

Type `T` can be a String or a Julia expression

## Examples
```julia
julia> mcall(\"integrate(sin(y)^2, y)\")

julia> mcall(:(integrate(1/(1+x^2), x)))

```
"""
mcall(expr::T) where T = mcall(MExpr(expr))

# Convert a Julia statement saved in `js` to a MExpr
mcallj(js::AbstractString) = mcall(MExpr(Meta.parse(js)))

function ==(m::MExpr, n::MExpr)
    r, s = normalize(m).strv("$x : $s")ec, normalize(n).strvec
    length(r) != length(s) && (return false)

    return all(j=>!occursin("false", mcall("is($(r[j]) = $(s[j]))")), 1:length(r))
end

macro mx(x, args...)
    r = MExpr[]
    if x == :use
        @assert !isempty(args) "Incomplete arguments for @mx"
        for a in args[1:end]
            s = jexpr_to_mstr(Main.eval(a))
            push!(r, mcall("$a : $s"))
        end
        return (length(r) == 1 ? r[1] : r)
    elseif x == :julia_export
        # TODO: not implemented yet
        @assert !isempty(args) "Incomplete arguments for @mx"
        x = args[1]
        return mcall(MExpr(x))
    end
    isempty(args) && return mcall(MExpr(x))

    # each element can be a Maxima statement, or a Maxima equation label
    argv = append!(Any[x], args)
    arglen = length(argv)
    i = 1
    while i <= arglen
        if argv[i].head === :(=)
            if i+1 <= arglen && argv[i+1] isa Symbol
                # equation definition
                push!(r, mcall(MExpr(string(argv[i+1])*": "*jexpr_to_mstr(x, defeqn=true))))
                i += 2
            else
                # variable definition
                push!(r, mcall(jexpr_to_mstr(x)))
                i += 1
            end
        else
            # function call or others
            push!(r, mcall(jexpr_to_mstr(x)))
            i += 1
        end
    end
    return (length(r) == 1 ? r[1] : r)
end

