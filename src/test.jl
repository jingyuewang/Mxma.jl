using MathMLRepr
using Mxma

#m = MExpr(:(f(x)=x^2+1))
#m = MExpr(:(a=1; b=x^2+2x+1))
#r = mcall1(m)
#println(r)

# Julia expressions => mstr
#mcall(:(b=x^2+2x+1))
#m1 = MExpr("sprint(b)\$")
#m1 = MExpr("sprint(b);")
#m1 = MExpr("a: 1\$ c: 2\$ print(b)\$")
#m1 = MExpr("sprint(b)")
#m1 = MExpr("a: matrix([1,2],[3,4])")
#mcall(MExpr(:(A = [1 2; 3 4])))
#mcall(MExpr(:(A = [1/2, 2/3, 3/4, 4/5])))
#mcall(MExpr(:(A = 34//87)))
#m1 = MExpr("sprint(A)\$")
#r = mcall(m1)
#println(r)

# Julia types => mstr
#A = rand(1:9, 2, 4)
#A = [1//2, 3//4, 5//6]
#A = rand(1:9, 3, 3); A = convert.(Rational, A)
#A = [:π :(x^2); :(x+1) :(2x)]

#s = @mxset(A)
# Some examples
#s = @mxc f = x^2+2x+1
#println(s)
#r = @mxc factor(f)
#println(r)
#r = @mxc expand((x+1)*(x+2))
#println(r)
#A = rand(1:9, 2, 4)
#r = @mxset A
#println(r)
#r = @mxc 2*[1,2,3,4]
#println(r)

# latex output
#r = @mxc f = x^2 + 3x + 2
#r = @mxc factor(f)
#io = IOBuffer()
#show(io, MIME("text/mathml+xml"), r)
#s = String(take!(io))
#println(s)

# Issue
#@mxc begin a=1;b=2 end
#r = @mx b=[1/2 2/3; 3/4 4/5]
#r = @mx b=[1 2; 3 4]
#r = @mx "cauchy_matrix([a,b,c])"
#io = IOBuffer()
#show(io, MIME("text/mathml+xml"), r)
#s = String(take!(io))
#println(s)
#println(typeof(r))
#println(length(Base.Multimedia.displays))
#display(r);println()
#display([r, r])
#display("text/plain", r); println()

# convert MExpr to LATEX
s = @mx fex=t^(1/2)*e^(-a*t/4)
t = convert(LATEX, s)
println(s)
println(t)
