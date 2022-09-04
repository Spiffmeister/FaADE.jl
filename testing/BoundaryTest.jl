


push!(LOAD_PATH,".")

using SBP_operators




Δx = 0.5
g₀(t) = 0.0




SAT_DL = Boundary_Dirichlet(g₀,Δx,Left,1,2)
SAT_DR = Boundary_Dirichlet(g₁,Δx,Left,1,2)



