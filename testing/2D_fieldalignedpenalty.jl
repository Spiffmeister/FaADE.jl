"""
    Testing using a parallel penalty where the field is aligned with the mesh, i.e. the grid point
    ``B:(x_i,y_j)\to (x_i,y_j)`` 
"""



cd("..")
push!(LOAD_PATH,"./SBP_operators")
push!(LOAD_PATH,"./plas_diff")

𝒟x = [0.0,1.0]
𝒟y = [0.0,1.0]

nx = ny = 41

Dom = Grid2D(𝒟x,𝒟y,nx,ny)

kx = ones(Float64,nx,ny)
ky = ones(Float64,nx,ny)

Δt = 1.0*min(Dom.Δx^2,Dom.Δy^2)
t_f = 100Δt

u₀(x,y) = x

BoundaryLeft = Boundary(Dirichlet,t->0.0,Left,1)
BoundaryRight = Boundary(Dirichlet,t->1.0,Right,1)
BoundaryUpDown = PeriodicBoundary(2)

order = 2
method = :cgie

P = VariableCoefficientPDE2D(u₀,kx,ky,BoundaryLeft,BoundaryRight,BoundaryUpDown)

H_x = SBP_operators.build_H(ny,order)
H_y = SBP_operators.build_H(nx,order)

H_x = 1.0 ./H_x.^2
H_y = 1.0 ./H_y.^2


κ_para = 1.0
τ_para = -1.0


function penalty_fn(u,u₀,Δt)
    for j = 1:ny
        for i = 1:nx
            u[i,j] = 1.0/(1.0 - κ_para/2.0 * Δt * (H_y[i] + H_x[j])) * (u₀[i,j] - Δt*τ_para/4.0 * (H_y[i]+H_x[j])*(u₀[i,j] + u₀[i,j]))
        end
    end
end
