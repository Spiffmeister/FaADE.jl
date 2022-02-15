module SBP_operators

    # Inbuild julia packages
    using LinearAlgebra


    # Include the files (add as completed)
    include("Op_Dx.jl")
    include("Op_Dxx.jl")


    # Export the functions for direct user interaction
    export Dₓ, Dₓₓ, Dₓ!, Dₓₓ!

end # module
