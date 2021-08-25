# 	This file is part of Maxima.jl. It is licensed under the MIT license
#   Copyright (c) 2016 Nathan Smith
#                 2021 Jingyue Wang

export string,
       show

import Base: string,
             show
#import Base.Multimedia: display

string(m::MExpr) = convert(String, m)

show(io::IO, m::MExpr) = print(io, convert(String, m))
function show(io::IO, ::MIME"text/plain", m::MExpr)
    check = "'(" * replace(convert(String, m), r";" => ");\n'(") * ")"
    write(ms.input, "$(replace(check, r";" => "\$"))\$\n print(ascii(4))\$")
    out = (readuntil(ms.output, EOT) |> String
                                     |> str -> rstrip(str, EOT))
    if occursin(synerr, out) || occursin(maxerr, out)
        @warn "Invalid Maxima expression"
        print(io, out)
    else
        mcall(mx"display2d: true")
        write(ms, replace(check, r";" => "; print(ascii(3))\$ "))
        str = rstrip(read(ms), '\n')
        mcall(mx"display2d: false")
        sp = split(str, ETX)
        for k in sp; print(io, k); end
    end
end


function show(io::IO, ::MIME"text/latex", m::MExpr)
    check = "'(" * replace(convert(String, m), r";" => ")\$\n'(") * ")"
    write(ms.input, "$check\$\n print(ascii(4))\$")
    out = (readuntil(ms.output, EOT) |> String
                                     |> str -> rstrip(str, EOT))
    if occursin(synerr, out) || occursin(maxerr, out)
        @warn "Invalid Maxima expression"
        print(io, out)
    else
        write(ms, "tex1('("*replace(convert(String, m), r";" => "))\$\ntex1('(") * "))")
        texstr = read(ms)
        print(io, replace(texstr, r"\nfalse\n" => ""))
    end
end

# A second approach for this function, which probably is better, is first
# converting `m` into a Julia expression, then use `MathMLRepr.jl`
function show(io::IO, ::MIME"text/mathml+xml", m::MExpr)
    check = "'(" * replace(convert(String, m), r";" => ")\$\n'(") * ")"
    write(ms.input, "$check\$\n print(ascii(4))\$")
    out = (readuntil(ms.output, EOT) |> String
                                     |> str -> rstrip(str, EOT))
    if occursin(synerr, out) || occursin(maxerr, out)
        @warn "Invalid Maxima expression"
        print(io, "Invalid Maxima Expression\n")
        print(io, out)
    else
        write(ms, "tex1('("*replace(convert(String, m), r";" => "))\$\ntex1('(") * "))")
        texstr = read(ms)
        s = replace(texstr, r"\nfalse\n" => "")
        # Warning! quick hacks inside this function
        s = clean_matrix_env(s)
        ## debug
        #write(io, s)
        ## Translate s into MathML.
        ts = TokenStream(s)
        nodes = exlist_to_mml(parse!(ts))
        if ts.top !== nothing && ts.top.val == "inline"
            o = mml_element(:math, nodes) 
        else
            o = mml_element_attr(:math, nodes, "display=\"block\"")
        end
        write(io, o)
    end
end


