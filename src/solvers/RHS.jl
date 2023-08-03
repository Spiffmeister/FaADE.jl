

struct Derivative{OrderX,OrderY,MixedOrderX,MixedOrderY} end



"""
    CGRHS!

"""
function CGRHS! end
function CGRHS!(cache::AT,CG::ConjGradBlock{TT,1,AT},D::newDataBlock) where {TT,AT}
    D₂!(cache,D.Data.u,D.Data.K,G.nx,G.ny,G.Δx,G.Δy,D.order,D.order)
    commBoundaries(DBlock,G,t,Δt)
    applySATs(CG,D,SolutionMode)
end
function CGRHS!(cache::AT,CG::ConjGradBlock{TT,2,AT},D::newDataBlock) where {TT,AT}
    D₂!(cache,D.Data.u,D.Data.K,G.nx,G.ny,G.Δx,G.Δy,D.order,D.order)
    commBoundaries(DBlock,G,t,Δt)
    applySATs(CG,D,SolutionMode)
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

