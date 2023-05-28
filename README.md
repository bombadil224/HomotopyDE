# HomotopyDE
Solve Differential equations using the Homotopy method.
Currently capable of basic ODE, PDE and systems of PDE’s. 

# Example
basic example of solving the wave equation

begin
    using ModelingToolkit
    using HomotopyDE

    println()
    println("wave equation")

    @variables t x ħ f(..)
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

    F = HomotopyDE.homotopy_solver(pde_system, degree, order, initialize)
   
    phi = HomotopyDE.power_series(F)
    phi = Symbolics.substitute(phi, ħ=>-1)

    println(phi)
end

[![Build Status](https://github.com/bombadil224/HomotopyDE.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/bombadil224/HomotopyDE.jl/actions/workflows/CI.yml?query=branch%3Amain)
