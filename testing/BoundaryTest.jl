


push!(LOAD_PATH,".")

using SBP_operators




g₀(t) = 0.0


SAT₀ = Boundary(Dirichlet,g₀,Left,order)


