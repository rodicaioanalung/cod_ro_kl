
using Plots
using LaTeXStrings


f(x) = exp(x - 2) - x

function g(x)
    return  exp(x - 2) - x
end

f(2)
g(2)


a, b = -2, 6

x = a:0.01:b;
y = f.(x);

plot(x, y, label=L"e^{x-2}-x")

# functions


""" 
An algorithm for bracketing an interval in which a local minimum must exist.
	 It takes as input a univariate function fand starting position x, 
		which defaults to 0. The starting step size s and the expansion factor
		 k can be specified. It returns a tuple containing the new interval [a, b]. 
"""
function bracket_minimum(f, x=0; s=1e-2, k=2.0)
    a, ya = x, f(x)
    b, yb = a + s, f(a + s)
    if yb > ya
        a, b = b, a
        ya, yb = yb, ya
        s = -s
    end
    while true
        c, yc = b + s, f(b + s)
        if yc > yb
            return a < c ? (a, c) : (c, a)
        end
        a, ya, b, yb = b, yb, c, yc
        s *= k

    end
end

print(bracket_minimum(f, 3; s=0.1, k=1.5))

using Base.MathConstants

"""
Fibonacci search
to be run on univariate func-
tion f, with bracketing interval
[a, b], for n > 1 function evalua-
tions. It returns the new interval
(a, b). The optional parameter ϵ
controls the lowest-level interval.
The golden ratio φ is defined in
Base.MathConstants.jl.
"""
function fibonacci_search(f, a, b, n; ϵ=0.00001)
    s = (1 - √5) / (1 + √5)
    ρ = 1 / (φ * (1 - s^(n + 1)) / (1 - s^n))
    r = (1 - ρ) * a + ρ * b # right sample point
    yr, n = f(r), n - 1
    while n > 0
        if n > 1
            l = ρ * a + (1 - ρ) * b
        else
            l = ϵ * a + (1 - ϵ) * r
        end
        ρ = 1 / (φ * (1 - s^(n + 1)) / (1 - s^n))
        yl, n = f(l), n - 1
        if yl < yr
            r, b, yr = l, r, yl
        else
            a, b = b, l
        end
        println("n=$(n), $(a), $(b)")
    end
    return a < b ? (a, b) : (b, a)

end

fibonacci_search(f, a, b, 5)


"""
Golden section
search to be run on a univariate
function f, with bracketing inter-
val [a, b], for n > 1 function eval-
uations. It returns the new inter-
val (a, b). Julia already has the
golden ratio φ defined. Guaran-
teeing convergence to within ϵ re-
quires n = (b−a)/(ϵ ln ϕ) itera-
tions.
"""
function golden_section_search(f, a, b, n)
    ρ = φ - 1
    d = ρ * b + (1 - ρ) * a
    yd = f(d)
    for i ∈ 1:(n-1)
        c = ρ * a + (1 - ρ) * b
        yc = f(c)
        if yc < yd
            b, d, yd = d, c, yc
        else
            a, b = b, c
        end
        println("i=$(i), $(a), $(b)")
    end
    return a < b ? (a, b) : (b, a)

end

golden_section_search(f, a,b,25)


"""
Quadratic fit search
to be run on univariate function f,
with bracketing interval [a, c] with
a< b< c. The method will run for
n function evaluations. It returns
the new bracketing values as a tu-
ple, (a, b, c).
"""
function quadratic_fit_search(f, a, b, c, n)
    ya, yb, yc = f(a), f(b), f(c)
    for i in 1:(n-3)
        x = 0.5 * (ya * (b^2 - c^2) + yb * (c^2 - a^2) + yc * (a^2 - b^2)) /
            (ya * (b - c) + yb * (c - a) + yc * (a - b))
        yx = f(x)
        if x > b
            if yx > yb
                c, yc = x, yx
            else
                a, ya, b, yb = b, yb, x, yx
            end
        elseif x < b
            if yx > yb
                a, ya = x, yx
            else
                c, yc, b, yb = b, yb, x, yx
            end

        end
    end
    return (a, b, c)
end

"""
	bisection(f′, a, b, ϵ)

 The bisection al-
gorithm, where f′ is the deriva-
tive of the univariate function we
seek to optimize. We have a < b
that bracket a zero of f′. The in-
terval width tolerance is ϵ. Calling
bisection returns the new brack-
eted interval [a, b] as a tuple.
The prime character′is not an
apostrophe. Thus, f′is a variable
name rather than a transposed vec-
tor f. The symbol can be created by
typing prime and hitting tab.
"""
function bisection(f′, a, b, ϵ)
    if a > b

        a, b = b, a
    end # ensure a < b
    ya, yb = f′(a), f′(b)
    if ya == 0

        b = a
    end
    if yb == 0

        a = b
    end
    while b - a > ϵ
        x = (a + b) / 2
        y = f′(x)
        if y == 0
            a, b = x, x
        elseif sign(y) == sign(ya)
            a = x
        else
            b = x
        end
    end
    return (a, b)
end

g(x) = x^3-x
x = 0.5:0.1:2
y = g.(x);
plot(x, y, label = "g")


bisection(g,0.5,2, 0.000001)


"""
An algorithm for
finding an interval in which a sign
change occurs. The inputs are the
real-valued function f′defined on
the real numbers, and starting in-
terval [a, b]. It returns the new in-
terval as a tuple by expanding the
interval width until there is a sign
change between the function eval-
uated at the interval bounds. The
expansion factor kdefaults to 2.
"""
function bracket_sign_change(f′, a, b; k=2)
    if a > b

        a, b = b, a
    end # ensure a < b
    center, half_width = (b + a) / 2, (b - a) / 2
    while f′(a) * f′(b) > 0
        half_width *= k
        a = center - half_width
        b = center + half_width
    end
    return (a, b)

end


bracket_sign_change(g, 0.7, 0.9)


# functions to plot : 

# f(x) = sin(x)
# f(x) = x^2
# f(x) = -x^2
# f(x,y) = x^2 + y^2
# f(x,y) = x^2 - y^2