using HomotopyDE
using Test
using Symbolics
using ModelingToolkit

function compair_F(test, Range, X)
err = []
    
for i in Range[1][1]:.1:Range[1][2]
    for j in Range[2][1]:.1:Range[2][2]

        test = Symbolics.substitute(test, X[1]=>i)
        test = Symbolics.substitute(test, X[2]=>j)

        push!(err, test)
        end
    end
return(maximum(err))
end

@testset "HomotopyDE.jl" begin
    # Write your tests here.
    @parameters x C
    P = 1 + x + 2 * x^2 + C*x^3
    F = HomotopyDE.Integrate_Polynomial(P, x)
    Dx = Symbolics.Differential(x)
    @test isequal(P, Symbolics.expand_derivatives(Dx(F)))



    begin
        println()
        println("Homogeneous Nonlinear system")

        @variables t x ħ u(..) v(..)
        Dt = Differential(t)
        Dx = Differential(x)

        ###Homogeneous Non-Linear system of PDEs
        eq  = [Dt(u(t, x)) + v(t,x) * (Dx(u(t, x))) + u(t, x) - 1 ~ 0,
              Dt(v(t, x)) + u(t,x) * (Dx(v(t, x))) - v(t, x) + 1 ~ 0 ]

        println(eq)

        degree = 1

       initialize = [[exp(x)], [exp(-x)]]

        @named NL_system = NonlinearSystem(eq, [u(t,x), v(t,x)], [t,x])

        order = 7
    end

    F = homotopy_solver(NL_system, degree, order, initialize)

    u = power_series(F[1])
    v = power_series(F[2])

    u = Symbolics.substitute(u, ħ=>-1)
    v = Symbolics.substitute(v, ħ=>-1)

    test = abs(exp(x-t) - u) / exp(x-t) * 100
    Range = [[2, 3], [2, 3]]

    @test compair_F(test, Range, NL_system.ps) ≈ 0 atol=4.0

    test = abs(exp(-x+t) - u) / exp(-x+t) * 100
    @test compair_F(test, Range, NL_system.ps) ≈ 0 atol=4.0
end