"""
Pass this the function that needs solved.
    The degree of the Differential Equation.
    The order of the solution that you want.
"""
function homotopy_solver(sys::ODESystem, degree::Int, order = 4::Int, initialize=Nothing)

    u = []      # this is the power series terms that we need to find
    R = []          # this is the Deformation equation terms
    A = []  
    phi = 0     # this is the series solution of the equation

    # now define the oporator A
    # from the differential equation A(f) = 0
    push!(A, sys.eqs[1].lhs - sys.eqs[1].rhs)
    push!(u, U_0(sys, degree, initialize))

    for i in 1:order
        R = Deformation_eqn(A, [u], sys.states)
        push!(u, Next_U(R, [u], t, degree)[1])

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


function U_0(sys::ODESystem, degree, initialize)
    @parameters C[1:degree]
    Result = 0
    
    if initialize == Nothing
        #Result = series_solution(sys.iv, degree)
        for i in 1:degree
            Result = Result + C[i] * sys.iv^(i-1)
        end
    else
        Result = initialize
        for i in 2:size(initialize)[1]
            Result = Result + (1 / (i-1)) * initialize[i] * t^(i - 1)
        end
    end
    return(Result)
end


function solve_boundary(u, sys, BC) # , degree)

    @parameters C[1:size(BC)[1]]
    Dt = Differential(sys.iv)
    boundary = []
    phi = 0
    for n in 1:size(u)[1]
        phi = phi + u[n]
    end

    for i in 1:size(BC)[1]
        push!(boundary, Symbolics.substitute(phi, t=>BC[i][1]) ~ BC[i][2])
        phi = Symbolics.expand_derivatives(Dt(phi))
    end

    Result = Symbolics.solve_for(boundary, C)

    for i in 1:size(u)[1]
        for j in 1:size(Result)[1]
            u[i] = Symbolics.substitute(u[i], C[j]=>Result[j])
        end
    end
    return u
end
