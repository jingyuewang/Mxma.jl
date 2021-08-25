
abstract type SMExpr end
abstract type Substitutable end

Numeric = Union{Integer, AbstractFloat, Nothing}

struct Num <: SMExpr
    literal::String
    num::Numeric
end

struct Var <: Substitutable
    id::String
end

struct FunExp <: Substitutable
end
