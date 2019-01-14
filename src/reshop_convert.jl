# Converts terms in the expression tree to supported expressions
convert_formula(c) = c
convert_formula(c::LinearityExpr) = convert_formula(c.c)

function convert_formula(c::Expr)
    map!(convert_formula, c.args, c.args)

    if c.head == :comparison
        n = length(c.args)
        @assert isodd(n)

        if VERSION >= v"0.5-"
            # All :comparison now have more than two comparisons.
            # When this becomes default, the `n > 3` test below can be removed
            @assert(n > 3)
        end

        # If more than binary comparison, we need to chain them together
        # Get binary comparisons in sequence and chain with nested &&
        # It looks like we might be able to use an n-ary &&, but that's not how
        # Julia parses chained && expressions.
        if n > 3
            new_expr = extract_binary_comparison(c, 1)
            for i in 3:2:(n - 2)
                new_expr = Expr(:&&, new_expr, extract_binary_comparison(c, i))
            end
            c = new_expr
        end

    elseif c.head == :call
        op = c.args[1]
        # TODO(xhub) we support n-ary multiplication
        # Distinguish unary, binary and n-ary plus
        if op == :+
            if length(c.args) > 3
                # n-ary plus to sum
                c.args[1] = :sum
            elseif length(c.args) == 2
                # Remove unary plus
                c = c.args[2]
            end
        # Distinguish unary and binary minus
        elseif op == :-
            n = length(c.args)
            if n == 2
                if c.args[2] == 0
                    c = 0
                else
                    c.args[1] = :neg
                end
            end
        # TODO(xhub) we support n-ary multiplication
        # .nl has no n-ary multiplication so we need to convert to binary
        elseif op == :*
            # Loop for each arg after the 3rd (more than binary)
            for _ in 1:(length(c.args) - 3)
                # Combine last term with previous to form a binary * expression
                arg = pop!(c.args)
                c.args[end] = :($(c.args[end]) * $arg)
            end
        # Handle normal conversion cases
        else
            if length(c.args) == 2
                c = get(unary_special_cases, c.args[1],
                        (x) -> Expr(:call, c.args[1], x))(c.args[2])
            end
        end
    end
    c
end

# Extracts `expression relation expression` from a larger comparison expression
function extract_binary_comparison(c::Expr, start::Integer)
    @assert c.head == :comparison
    @assert start <= length(c.args) - 2

    if VERSION < v"0.5-"
        Expr(:comparison, c.args[start], c.args[start + 1], c.args[start + 2])
    else
        Expr(:call, c.args[start + 1], c.args[start], c.args[start + 2])
    end
end

const unary_special_cases = Dict(
:cbrt  => (x) -> :($x ^ (1 / 3)),
:abs2  => (x) -> :($x ^ 2),
:inv   => (x) -> :(1 / $x),
:log2  => (x) -> :(log($x) / log(2)),
:log1p => (x) -> :(log(1 + $x)),
:exp2  => (x) -> :(2 ^ $x),
:expm1 => (x) -> :(exp($x) - 1),

:sec   => (x) -> :(1 / cos($x)),
:csc   => (x) -> :(1 / sin($x)),
:cot   => (x) -> :(1 / tan($x)),

:asec  => (x) -> :(acos(1 / $x)),
:acsc  => (x) -> :(asin(1 / $x)),
:acot  => (x) -> :(pi / 2 - atan($x)),

:sind  => (x) -> :(sin(pi /180 * $x)),
:cosd  => (x) -> :(cos(pi /180 * $x)),
:tand  => (x) -> :(tan(pi /180 * $x)),
:secd  => (x) -> :(1 / cos(pi / 180 * $x)),
:cscd  => (x) -> :(1 / sin(pi / 180 * $x)),
:cotd  => (x) -> :(1 / tan(pi / 180 * $x)),

:asind => (x) -> :(asin($x) * 180 / pi),
:acosd => (x) -> :(acos($x) * 180 / pi),
:atand => (x) -> :(atan($x) * 180 / pi),
:asecd => (x) -> :(acos(1 / $x) * 180 / pi),
:acscd => (x) -> :(asin(1 / $x) * 180 / pi),
:acotd => (x) -> :((pi / 2 - atan($x)) * 180 / pi),

:sech  => (x) -> :(1 / cosh($x)),
:csch  => (x) -> :(1 / sinh($x)),
:coth  => (x) -> :(1 / tanh($x)),

:asech => (x) -> :(acosh(1 / $x)),
:acsch => (x) -> :(asinh(1 / $x)),
:acoth => (x) -> :(atanh(1 / $x)),
)
