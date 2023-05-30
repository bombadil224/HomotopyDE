using HomotopyDE
using Test
using Symbolics
using ModelingToolkit

@testset "HomotopyDE.jl" begin
    # Write your tests here.
    @parameters x C
    P = 1 + x + 2 * x^2 + C*x^3
    F = HomotopyDE.Integrate_Polynomial(P, x)
    Dx = Symbolics.Differential(x)
    @test isequal(P, Symbolics.expand_derivatives(Dx(F)))


end
