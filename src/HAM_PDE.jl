function homotopy_solver(sys::PDESystem, degree::Int, order = 4::Int, initialize=Nothing)

    u = []      # this is the power series terms that we need to find
    R = []          # this is the Deformation equation terms
    phi = 0     # this is the series solution of the equation
    A = []

    # now define the oporator A
    # from the differential equation A(f) = 0
    push!(A, sys.eqs[1].lhs - sys.eqs[1].rhs)
    push!(u, U_0(sys.ivs, degree, initialize))

    for i in 1:order
        R = Deformation_eqn(A, [u], sys.dvs)
        push!(u, Next_U(R, [u], sys.ivs[1], degree)[1])
  
        print("\r")
        print("Solved order ")
        print(i)
        print(" of ")
        print(order)
    end
    println("")
    println("done")
    return(u)
end