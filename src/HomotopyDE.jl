module HomotopyDE
    using Symbolics
    using ModelingToolkit


    include("HAM.jl")
    include("HAM_ODE.jl")
    include("HAM_PDE.jl")
    include("HAM_NLSys.jl")

    export homotopy_solver
    export power_series
    export solve_boundary

end
