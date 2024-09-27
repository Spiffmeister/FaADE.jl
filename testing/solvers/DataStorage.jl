


Dl = FaADE.SATs.SAT_Dirichlet(BxL,Dom1V.Δx,Left,1,order)
Dr = FaADE.SATs.SAT_Dirichlet(BxR,Dom1V.Δx,Right,1,order)
BD1V = (Dl,Dr)
P1V = newProblem1D(order,u₀,K,Dom1V,BD1V,F,nothing)


FaADE.solvers.DataMultiBlock(P::FaADE.Problems.Problem1D,G::FaADE.Helpers.GridMultiBlock,Δt::Float64,t::Float64) = FaADE.solvers.DataMultiBlock(P,G,Δt,t,0.0)

##======##
# newBoundaryData
##======##


Dom = Grid1D([0.0,1.0],101)

order = 2

K = 1.0

u₀(x) = 0.0

Dl = FaADE.SATs.SAT_Dirichlet(t->0,Dom.Δx,Left,1,order)
Dr = FaADE.SATs.SAT_Dirichlet(t->0,Dom.Δx,Right,1,order)
BD = (Dl,Dr)
P = newProblem1D(order,u₀,K,Dom,BD,nothing,nothing)

SC = FaADE.solvers.StepConfig{Float64}()

FaADE.solvers.newLocalDataBlock(P,Dom,SC)




##======##
# newInterfaceBoundaryData
##======##



##======##
# newLocalDataBlock
##======##


