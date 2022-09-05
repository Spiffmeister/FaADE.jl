


push!(LOAD_PATH,".")

using SBP_operators




Δx = 0.5
g₀(t) = 0.0



Gx = Grid1D([0.0,1.0], 41, 2)
Gy = Grid1D([0.5,0.68], 41, 2)

# SquareGrid



SAT_DL = Boundary_Dirichlet(g₀,Δx,Left,1,2)
SAT_DR = Boundary_Dirichlet(g₁,Δx,Right,1,2)



