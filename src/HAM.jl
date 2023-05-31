@parameters ħ q
Dq = Differential(q)
@variables t f(..)

function series_solution(X, degree)

    A = zeros(Int32, size(X)[1])
    answer = []
    println(A)

    while A[size(A)[1]] <= degree
        if sum(A) <= degree
            push!(answer, copy(A))
        end

        A[1] = A[1] + 1

        for i in 1:(size(A)[1] - 1)
            if sum(A) > degree
                A[i] = 0
                A[i+1] = A[i+1] + 1
            end
        end
    end

    @parameters C[1:(size(answer)[1])]

    result = 1
    total = 0
    for j in 1:size(answer)[1]
        for i in 1:size(pde_system.ivs)[1]
            result = result * X[i]^(answer[j][i])
        end
        total = total + C[j] * result
        result = 1
    end
    return(total)
end


function U_0(X, degree, initialize)
    @parameters C[1:degree]
    Result = 0
    
    if initialize == Nothing
        Result = series_solution(X, degree)
    else 
        Result = initialize[1]
        for i in 2:size(initialize)[1]
            Result = Result + (1 / (i-1)) * initialize[i] * t^(i - 1)
        end
    end

    return(Result)
end


"""
This function repetadly applys integration to
    the deformation equation to find the next
    term in the homotopy exspantion
"""
function L_inv(Q, x, degree, debug = 0)
    Result = Q
    for i in 1:degree
       # Result = integrate(Result, x)[1]
       Result = Integrate_Polynomial(Result, x, debug)
    end
    return(Result)
end


"""
This is a vary primetive integrator and needs imporuvment
    as it is the limiting factor in the functions that can be
    used in the equation being solved.
This function needs improvment
"""
function Integrate_Polynomial(Q, x, debug = 0)
    # at the vary least this function needs error checking
    Q = expand_derivatives(Q)
    Q = Symbolics.simplify(Q)
    Q = Symbolics.expand(Q)
    Result = 0

    if debug == 1
        render(latexify(Q))
    end

    while(Symbolics.degree(Q, x) != 0)
        n = Symbolics.degree(Q, x)
        Result = Result + Symbolics.coeff(Q, x^n) * x^(n+1) / (n+1)
        Q = Q - Symbolics.coeff(Q, x^n) * x^n

        if debug == 1
            println("the valu of n is")
            println(n)
            println("the valu of Q is")
            println(Q)
            println("the coeficent is")
            println(Symbolics.coeff(Q, x^n))
            println("the result is")
            println(Result)
        end
        Q = Symbolics.simplify(Q)
        Q = Symbolics.expand(Q)
    end

    Result = Result + Q * x
    return(Result)
end


"""
This finds the deformation equation of the oporator A
"""
function Deformation_eqn(A, u, X)
    phi = 0
    Result = []

    A_2 = deepcopy(A)

    for n in 1:size(u)[1]
        phi = power_series(u[n], q)
        for i in 1:size(A_2)[1]
            A_2[i] = Symbolics.substitute(A_2[i], X[n]=>phi)
        end
    end
    for i in 1:size(A_2)[1]
        for n in 1:size(u[i])[1] - 1
            A_2[i] = Dq(A_2[i]) / n
        end

        A_2[i] = Symbolics.expand_derivatives(A_2[i])
        A_2[i] = Symbolics.substitute(A_2[i], q=>0)
    end
    return(A_2)
end


"""
equation 18
solve for the next u_m needed for the series solution
    Apply the Invers Linear opporator to the deformation
    equation (integration)
"""
function Next_U(R, u, t, degree)
        Result = []

         # if this is not the first term of u then the last term
        # of u must be added to the term for the next u
        for i in 1:size(u)[1]
            if size(u[i])[1] > 1   # apply the Invers Linear opporator to the deformation equation
                push!(Result, last(u[i]) + ħ * L_inv(R[i], t, degree))
            else
                push!(Result, ħ * L_inv(R[i], t, degree))
            end
            Symbolics.simplify(Result[i])
            Symbolics.expand(Result[i])
        end
        return(Result)
    end


"""
This funciton makes a power series
    out of the terms that it is given in powers of x
"""
function power_series(C, x=1)
    phi = 0
    for i in 1:size(C)[1]
        phi = phi + C[i] * x^(i-1)
    end

    return(phi)
end
