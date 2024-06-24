using LinearAlgebra
using Revise
using FaADE


plot = true

θ = 0.5
order = 4

k_para = 1.0e6
k_perp = 1.0



# Domain
𝒟x = [-0.5,0.5]
𝒟y = [-0.5,0.5]


nx = 11
ny = 21
D1 = Grid2D([-0.5,0.0],𝒟y,nx,ny)
D2 = Grid2D([0.0,0.5],𝒟y,nx,ny)

Dom = GridMultiBlock((D1,D2), 
        ((Joint(2,Right),),
        (Joint(1,Left),)))







# Magnetic field
Ψ(x,y) = cos(π*x)*cos(π*y)


function B(X,x,p,t)
    bn = π * sqrt(abs(cos(x[1]*π)*sin(x[2]*π))^2 + abs(sin(x[1]*π)*cos(x[2]*π))^2)
    X[1] = π*cos(π*x[1])*sin(π*x[2])/bn
    X[2] = -π*sin(π*x[1])*cos(π*x[2])/bn
    if (x[1] == 0.5) && (x[2] == 0.5)
        X[1] = 0.0
        X[2] = -1.0
    elseif (x[1] == 0.5) && (x[2] == -0.5)
        X[1] = -1.0
        X[2] = 0.0
    elseif (x[1] == -0.5) && (x[2] == -0.5)
        X[1] = 0.0
        X[2] = 1.0
    elseif (x[1] == -0.5) && (x[2] == 0.5)
        X[1] = 1.0
        X[2] = 0.0
    elseif (x[1] == 0.0) && (x[2] == 0.0)
        X[1] = 0.0
        X[2] = 0.0
    end
    # X[3] = 0.0
end
MagField(X,t) = [
    π*cos(π*X[1])*sin(π*X[2]),
    -π*sin(π*X[1])*cos(π*X[2]),
    0.0
]

# Exact solution
T(X,t) = (1.0 - exp(-2.0*π^2*k_perp*t) )*Ψ(X[1],X[2])/k_perp # k_perp = 1


# Initial condition
u₀(x,y) = T([x,y],0.0)
# u₀(x,y) = T([x,y],Inf)
# Source term
F(X,t) = 2π^2*cos(π*X[1])*cos(π*X[2])


coord = :Cartesian

# Homogeneous boundary conditions
Boundary1Left    = SAT_Dirichlet((y,t) -> 0.0, D1.Δx , Left,  order)
Boundary2Right   = SAT_Dirichlet((y,t) -> 0.0, D2.Δx , Right, order)

Boundary1Up      = SAT_Dirichlet((x,t) -> 0.0, D1.Δy , Up,    order)
Boundary1Down    = SAT_Dirichlet((x,t) -> 0.0, D1.Δy , Down,  order)

Boundary2Up      = SAT_Dirichlet((x,t) -> 0.0, D2.Δy , Up,    order)
Boundary2Down    = SAT_Dirichlet((x,t) -> 0.0, D2.Δy , Down,  order)

BC = Dict(1 => (Boundary1Left,Boundary1Up,Boundary1Down), 2 => (Boundary2Right,Boundary2Up,Boundary2Down))



# Time setup
Δt = 0.1D1.Δx^2
t_f = 0.1
nf = round(t_f/Δt)
Δt = t_f/nf




gdata   = construct_grid(B,Dom,[-1.0,1.0],ymode=:ignore)
gdata_chk = construct_grid(B,Dom,[-1.0,1.0],ymode=:ignore,interpmode=:bilinear)


PData   = FaADE.ParallelOperator.ParallelMultiBlock(gdata,Dom,order,κ=k_para,interpolant=:bicubic)




# Build PDE problem
P       = Problem2D(order,u₀,k_perp,k_perp,Dom,BC,F,PData)



soln = solve(P,Dom,Δt,1.1Δt,solver=:theta,  θ=θ)
# soln = solve(P,Dom,Δt,t_f,  solver=:theta,  θ=θ)





###=== CHECKING ===###
D12 = Grid2D([-0.5,0.5],[-0.5,0.5],2nx-1,ny)
gdata_nn = construct_grid(B,D12,[-1.0,1.0],ymode=:ignore,interpmode=:nearest)
PData_nn = ParallelData(gdata_nn,D12,order,κ=k_para,interpolant=:nearest)

BC_nn = (
    SAT_Dirichlet((y,t) -> cos(0.5π)*cos(π*y[2]) , D12.Δx , Left,  order),
    SAT_Dirichlet((y,t) -> cos(-0.5π)*cos(π*y[2]), D12.Δx , Right, order),
    SAT_Dirichlet((x,t) -> cos(π*x[1])*cos(0.5π) , D12.Δy , Up,    order),
    SAT_Dirichlet((x,t) -> cos(π*x[1])*cos(-0.5π), D12.Δy , Down,  order)
)

P_nn = Problem2D(order,u₀,k_perp,k_perp,D12,BC_nn,F,PData_nn)

soln_nn = solve(P_nn,D12,Δt,1.1Δt,solver=:theta,θ=θ)




###=== CHECKING CHECKING CHECKING ===###
for i in 1:nx
    for j in 1:ny
        if Dom.Grids[gdata[1].Fplane.subgrid[i,j]][gdata[1].Fplane.x[i,j]] != D12[gdata_nn.Fplane.x[i,j]]
            error("Mismatch at ($(i),$(j)) in 1")
        end
    end
end
for i in 1:nx
    for j in 1:ny
        if Dom.Grids[gdata[2].Fplane.subgrid[i,j]][gdata[2].Fplane.x[i,j]] != D12[gdata_nn.Fplane.x[nx+i-1,j]]
            error("Mismatch at ($(i),$(j)) in 2)")
        end
    end
end
for i in 1:nx
    for j in 1:ny
        if Dom.Grids[gdata[1].Bplane.subgrid[i,j]][gdata[1].Bplane.x[i,j]] != D12[gdata_nn.Bplane.x[i,j]]
            error("Mismatch at ($(i),$(j)) in 1")
        end
    end
end
for i in 1:nx
    for j in 1:ny
        if Dom.Grids[gdata[2].Bplane.subgrid[i,j]][gdata[2].Bplane.x[i,j]] != D12[gdata_nn.Bplane.x[nx+i-1,j]]
            error("Mismatch at ($(i),$(j)) in 2)")
        end
    end
end




T_exact = [zeros(eltype(Dom.Grids[I]),size(Dom.Grids[I])) for I in eachindex(Dom.Grids)];
for I in eachindex(Dom.Grids)
    for J in eachindex(Dom.Grids[I])
        T_exact[I][J] = T(Dom.Grids[I][J],soln.t[2])
    end
end


# pollution = abs(1 - soln.u[2][floor(Int,nx/2)+1,floor(Int,ny/2)+1])
# rel_error = norm(T_exact .- soln.u[2])/norm(T_exact)
# println("n=",nx," poll=",pollution," relerr=",rel_error," abserr=",norm(T_exact .- soln.u[2])," t=",soln.t[2])
# @show norm(T_exact .- soln.u[2]) / sqrt(nx*ny)


if plot
    using GLMakie
    f = Figure(); 
    
    ax = Axis3(f[1,1]); 
    surface!(ax,Dom.Grids[1].gridx,Dom.Grids[1].gridy,T_exact[1])
    surface!(ax,Dom.Grids[2].gridx,Dom.Grids[2].gridy,T_exact[2])
    
    ax2 = Axis3(f[1,2]); 
    surface!(ax2,Dom.Grids[1].gridx,Dom.Grids[1].gridy,soln.u[2][1])
    surface!(ax2,Dom.Grids[2].gridx,Dom.Grids[2].gridy,soln.u[2][2])

    ax3 = Axis3(f[1,3]);
    surface!(ax3,D12.gridx,D12.gridy,soln_nn.u[2])


    # surface!(ax,Dom.gridx,Dom.gridy,soln.u[2])
    # surface!(ax,Dom.gridx,Dom.gridy,soln.u[2])
    # surface!(ax,Dom.gridx,Dom.gridy,w_f)
end


# f = Figure();
# ax = Axis(f[1,1])
# scatter!(ax,Dom.gridx[:],Dom.gridy[:])
# scatter!(ax,gdata.Fplane.x[:],gdata.Fplane.y[:],color=:red)


