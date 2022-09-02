


push!(LOAD_PATH,".")

using SBP_operators




g₀(t) = 0.0


SAT₀ = Boundary(Dirichlet,g₀,Left,order,Δx,dim=1)

SAT₁ = Boundary(Dirichlet,g₁,Right,order,Δx,dim=1)








