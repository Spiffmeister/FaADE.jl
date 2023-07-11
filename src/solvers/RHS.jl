

struct Derivative{OrderX,OrderY,MixedOrderX,MixedOrderY} end


struct RHS
    Diff    :: Derivative
    SATL    :: Union{SimultaneousApproximationTerm,nothing}
    SATR    :: Union{SimultaneousApproximationTerm,nothing}
    SATU    :: Union{SimultaneousApproximationTerm,nothing}
    SATD    :: Union{SimultaneousApproximationTerm,nothing}

    function RHS
    end
end




"""
    CGRHS!

"""
function CGRHS! end

function CGRHS!(::Derivative{2,2,0,0},RHS)
    D₂!(cache,u,Kx,Ky,nx,ny,Δx,Δy,orderx,ordery)

    ApplySAT(cache,u,Kx,RHS.SATL)
    ApplySAT(cache,u,Kx,RHS.SATR)
    ApplySAT(cache,u,Ky,RHS.SATU)
    ApplySAT(cache,u,Ky,RHS.SATD)
end



function CGRHS!(::Derivative{2,2,1,1})
end


"""
    ApplyBoundary
Apply the SATs or do nothing
"""
function ApplyBoundary end
"""
    ApplyBoundary(cache,u,k,SAT::SimultaneousApproximationTerm)
Apply the SAT to the right hand side
"""
function ApplyBoundary(SAT::SimultaneousApproximationTerm,cache,u,k)
    SAT(cache,u,k,SolutionMode)
end
"""
    ApplyBoundary
Apply the Periodic SAT to the right hand side
"""
function ApplyBoundary(SAT::SAT_Periodic,cache,u,k)
    SAT(cache,u,k)
end
"""
    ApplyBoundary(cache,u,k,SAT::nothing)
If there is no SAT, do nothing
"""
function ApplyBoundary(SAT::nothing,cache,u,k) end

