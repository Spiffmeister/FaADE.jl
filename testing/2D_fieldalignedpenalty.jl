"""
    Testing using a parallel penalty where the field is aligned with the mesh, i.e. the grid point
    ``B:(x_i,y_j)\to (x_i,y_j)`` 
"""



cd("..")
push!(LOAD_PATH,"./FaADE")
using FaADE

𝒟x = [0.0,1.0]
𝒟y = [0.0,1.0]

nx = ny = 41

Dom = Grid2D(𝒟x,𝒟y,nx,ny)

kx = ones(Float64,nx,ny)
ky = ones(Float64,nx,ny)

Δt = 1.0*min(Dom.Δx^2,Dom.Δy^2)
t_f = 1_000Δt

u₀(x,y) = x

BoundaryLeft = Boundary(Dirichlet,(t,y)->0.0,Left,1)
BoundaryRight = Boundary(Dirichlet,(t,y)->1.0,Right,1)
BoundaryUpDown = PeriodicBoundary(2)

order = 2
method = :cgie

P = VariableCoefficientPDE2D(u₀,kx,ky,order,BoundaryLeft,BoundaryRight,BoundaryUpDown)

H_x = FaADE.build_H(ny,order)
H_y = FaADE.build_H(nx,order)

H_x = 1.0 ./H_x.^2
H_y = 1.0 ./H_y.^2


κ_para = 1.0
τ_para = -1.0


function penalty_fn(u,u₀,Δt)
    for j = 1:ny
        for i = 1:nx
            u[i,j] = 1.0/(1.0 - κ_para*τ_para/2.0 * Δt * (H_y[i] + H_x[j])) * (u₀[i,j] - Δt*κ_para*τ_para/4.0 * (H_y[i]+H_x[j])*(u₀[i,j] + u₀[i,j]))
        end
    end
end


soln = solve(P,Dom,Δt,5.1Δt,:cgie,adaptive=true,penalty_func=penalty_fn)
soln = solve(P,Dom,Δt,t_f,:cgie,adaptive=true,penalty_func=penalty_fn)

println("Plotting")
using Plots
surface(soln.grid.gridy,soln.grid.gridx,soln.u[2],
    xlabel="y",ylabel="x",zlabel="Temp")