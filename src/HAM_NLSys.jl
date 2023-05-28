function homotopy_solver(sys::NonlinearSystem, degree::Int, order = 4::Int, initialize=Nothing)
    u = []      # this is the power series terms that we need to find
    R = []          # this is the Deformation equation terms
    #phi = []     # this is the series solution of the equation
    phi = 0     # this is the series solution of the equation

    for i in 1:size(sys.eqs)[1]
        push!(u,[])
    end

    A = []
    # now define the oporator A
    # from the differential equation A(f) = 0
    for i in 1:size(sys.eqs)[1]
        push!(A, sys.eqs[i].lhs - sys.eqs[i].rhs)
        push!(u[i], U_0(sys.ps[1], degree, initialize[i]))
    end


        for i in 1:order
            R = Deformation_eqn(A, u, sys.states)
            F = Next_U(R, u, sys.ps[1], degree)

            for j in 1:size(u)[1]
                push!(u[j], F[j])
            end

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
