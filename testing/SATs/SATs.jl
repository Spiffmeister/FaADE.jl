

using Test
using FaADE



𝒟 = [0.0,1.0]
n = 21
Dom = Grid1D(𝒟, n)

order = 2
K = ones(Float64,n)



##======##
# Dirichlet SATs
##======##

SAT_Dirichlet_Left  = FaADE.SATs.SAT_Dirichlet((x,t)->0.0,Dom.Δx,Left,  order)
SAT_Dirichlet_Right = FaADE.SATs.SAT_Dirichlet((x,t)->0.0,Dom.Δx,Right, order)
SAT_Dirichlet_Up    = FaADE.SATs.SAT_Dirichlet((x,t)->0.0,Dom.Δx,Up,    order)
SAT_Dirichlet_Down  = FaADE.SATs.SAT_Dirichlet((x,t)->0.0,Dom.Δx,Down,  order)







##======##
# Neumann SATs
##======##

SAT_Neumann_Left  = FaADE.SATs.SAT_Neumann((x,t)->0.0,Dom.Δx,Left,  order)
SAT_Neumann_Right = FaADE.SATs.SAT_Neumann((x,t)->0.0,Dom.Δx,Right, order)
SAT_Neumann_Up    = FaADE.SATs.SAT_Neumann((x,t)->0.0,Dom.Δx,Up,    order)
SAT_Neumann_Down  = FaADE.SATs.SAT_Neumann((x,t)->0.0,Dom.Δx,Down,  order)



##======##
# Robin SATs
##======##




##======##
# Periodic SATs
##======##





##======##
# Interface SATs
##======##





