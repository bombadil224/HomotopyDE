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


    println()
    println("nonlinear second-order IVP")
    @variables t f(t) ħ
    Dt = Differential(t)
    Dtt = Differential(t)^2

    eqs = [ Dtt(f) ~ - Dt(f)^2 ]
    println(eqs)
    degree = 2

    @named sys = ODESystem(eqs, t)

    u = homotopy_solver(sys, degree)

    BC = []
    push!(BC, [0,1])
    push!(BC, [0,2])

    u = solve_boundary(u, sys, BC)

    test = []
    push!(test, 1 + 2 * t)
    push!(test, 2 * t^2 * ħ)
    push!(test, 2 * t^2 * ħ + 2 * t^2 * ħ^2 + 8 * t^3 / 3 * ħ^2)
    push!(test, 2 * t^2 * ħ + 4 * t^2 * ħ^2 + 2 * ħ^3 * t^2 + 16 * t^3 / 3 * ħ^2
                + 16 * ħ^3 / 3 * t^3 + 4 * ħ^3 * t^4)

    test_score = 0
    for i in 1:size(test)[1]
        test_u = u[i] - test[i]
        test_u = Symbolics.expand(test_u)
        test_u = Symbolics.simplify(test_u)
        if isequal(test_u, 0) == 1
            test_score = test_score + 1
        else
            print("Test failed on eqution ")
            println(i)
        end
    end

    @test test_score == 4


    begin
    println()
    println("wave equation")

    @variables t x f(..) ħ
    Dxx = Differential(x)^2
    Dtt = Differential(t)^2
    Dt = Differential(t)^2

    ###2D PDE
    eq  = Dtt(f(t,x)) ~ Dxx(f(t,x))
    println(eq)

    degree = 2

    ### Initial and boundary conditions
    bcs = []

    ### Space and time domains
    domains = [t ∈ (0.0,1.0),
               x ∈ (0.0,1.0)]

    @named pde_system = PDESystem(eq,bcs,domains,[t,x],[f(t,x)])
    initialize = [x^2, 1]
    order = 3
end

    F = homotopy_solver(pde_system, degree, order, initialize)

    phi = power_series(F)
    phi = Symbolics.substitute(phi, ħ=>-1)

    @test isequal(phi, t + x^2 + t^2)


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