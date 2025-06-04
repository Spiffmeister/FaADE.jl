using LinearAlgebra
using Revise
using FaADE


plot = false
poincare = false
savefigs = false
save_reference_solution = false

order = 2


# Time setup
Δt = 1e-4
# Δt = 0.1
t_f = 1e-3

k_para = 1.0e8
k_perp = 1.0

#=== DOMAIN SETUP ===#
𝒟x = [0.0, 1.0]
𝒟y = [0.0, 2π]
nx = 201
ny = 41

Dx(x, nx) = sinh(0.15 * x * (nx / 51)^1.3) / 2sinh(0.15 * (nx / 51)^1.3) + 0.5
Dy(y) = y

𝒟x, 𝒟y = FaADE.Grid.meshgrid(Dx.(LinRange(-1.0, 1.0, nx), nx), Dy.(LinRange(0, 2π, ny)))
TestDom = Grid2D(𝒟x, 𝒟y, ymap=false)

D1 = Grid2D([0.0, 0.3], [0.0, 2π], 41, ny)
D2 = Grid2D([0.3, 0.65], [0.0, 2π], 201, ny)
D3 = Grid2D([0.65, 1.0], [0.0, 2π], 41, ny)

joints = ((Joint(2, Right),),
    (Joint(1, Left), Joint(3, Right)),
    (Joint(2, Left),))

Dom = GridMultiBlock((D1, D2, D3), joints)

u₀(x, y) = x

coord = :Cartesian


#=== BOUNDARY CONDITIONS ===#
Bl1 = SAT_Dirichlet((y, t) -> 0.0, D1.Δx, Left, order, D1.Δy, coord) #x=0
Br3 = SAT_Dirichlet((y, t) -> 1.0, D3.Δx, Right, order, D3.Δy, coord) #x=1

Bu1 = SAT_Periodic(D1.Δy, order, Up, D1.Δx, coord) #Block 1
Bd1 = SAT_Periodic(D1.Δy, order, Down, D1.Δx, coord) #Block 1

Bu2 = SAT_Periodic(D2.Δy, order, Up, D2.Δx, coord) #Block 2
Bd2 = SAT_Periodic(D2.Δy, order, Down, D2.Δx, coord) #Block 2

Bu3 = SAT_Periodic(D3.Δy, order, Up, D3.Δx, coord) #Block 3
Bd3 = SAT_Periodic(D3.Δy, order, Down, D3.Δx, coord) #Block 3

BC = Dict(1 => (Bl1, Bu1, Bd1), 2 => (Bu2, Bd2), 3 => (Br3, Bu3, Bd3))



#=== PARALLEL MAP===#
δ = 0.015
xₛ = 0.5

function B(X, x, p, t)
    # bn = 1 + abs( δ*x[1]*(x[1]-1)*sin(x[2]) )^2 + abs( 2*x[1] - 2*xₛ + δ*(1-x[1])*cos(x[2]) - δ*x[1]*cos(x[2]) )^2
    # bn = sqrt(bn)
    X[1] = δ * x[1] * (1 - x[1]) * sin(x[2])#/bn
    X[2] = (2x[1] - 2 * xₛ + δ * (1 - x[1]) * cos(x[2]) - δ * x[1] * cos(x[2]))#/bn
end


if poincare
    using GLMakie
    include("../../BADESBP_examples/FieldLines.jl")
    # include("../paper_JCP2023/FieldLines.jl")
    poindata = FieldLines.construct_poincare(B, [0.0, 1.0], [0, 2π])
    # h = Figure(); ax_h = Axis(h[1,1]);
    h = Figure(size=(1600, 1200), fontsize=20)
    ax_h = Axis(h[1, 1], ylabel="x", xlabel="y", ylabelsize=30, xlabelsize=30;)
    xlims!(ax_h, 0, 2π)
    ylims!(ax_h, 0, 1)
    scatter!(ax_h, poindata.θ[:], poindata.ψ[:], markersize=3.0, color=:black, alpha=0.4)#,xlims=(0,2π),ylims=(0,1))
    ax_h.yticks = 0.0:0.1:1.0

    # contour!(ax_h,Dom.gridy[1,:],Dom.gridx[:,1],soln.u[2]',levels=25,linewidth=3)
end


gridoptions = Dict("xbound" => [0.0, 1.0], "ybound" => [0.0, 2π], "remapping" => :idw)

gdata = construct_grid(B, Dom, [-2.0π, 2.0π], gridoptions=gridoptions)
# gdata   = construct_grid(B,Dom,[-2.0π,2.0π],interpmode=:bicubic)
PData = FaADE.ParallelOperator.ParallelMultiBlock(gdata, Dom, order, κ=k_para)


# S(X,t) = (1-X[1]^2)^2
# S(X,t) = (1-(X[1]-1)^2)^2
S = nothing

# Build PDE problem
P = Problem2D(order, u₀, k_perp, k_perp, Dom, BC, S, PData)

println("Solving")

soln = solve(P, Dom, Δt, 1.1Δt, solver=:theta, θ=θ)
soln = solve(P, Dom, Δt, t_f, solver=:theta, θ=θ)




if plot
    println("plotting")
    using GLMakie

    colourrange = (minimum(minimum.(soln.u[2])), maximum(maximum.(soln.u[2])))

    # using CairoMakie
    f = Figure()
    ax_f = Axis3(f[1, 1])

    surface!(ax_f, Dom.Grids[1].gridx, Dom.Grids[1].gridy, soln.u[2][1], colorrange=colourrange)
    surface!(ax_f, Dom.Grids[2].gridx, Dom.Grids[2].gridy, soln.u[2][2], colorrange=colourrange)
    surface!(ax_f, Dom.Grids[3].gridx, Dom.Grids[3].gridy, soln.u[2][3], colorrange=colourrange)



    gridfig = Figure()
    ax_g = Axis(gridfig[1, 1])
    scatter!(ax_g, Dom.Grids[1].gridx[:], Dom.Grids[1].gridy[:], markersize=3.0)
    scatter!(ax_g, Dom.Grids[2].gridx[:], Dom.Grids[2].gridy[:], markersize=3.0)
    scatter!(ax_g, Dom.Grids[3].gridx[:], Dom.Grids[3].gridy[:], markersize=3.0)


    # lines!(ax2_f,Dom.gridx[:,1],soln.u[2][:,floor(Int,ny/2)+1])

    # wireframe!(ax,Dom.gridx,Dom.gridy,soln.u[2])
    # contour3d!(ax,soln.u[2],levels=100)

    # g = Figure();
    # ax_g = Axis(g[1,1]);
    # contour3d!(ax_g,Dom.gridy[1,:],Dom.gridx[:,1],soln.u[2]',levels=100)
    # Colorbar(g[1,2])
end




#=
if savefigs
    f = Figure(); ax = Axis(f[1,1]); lines!(ax,Dom.gridx[:,1],soln.u[2][:,floor(Int,ny/2)+1])
    save("SingleIsland_out/SingleIsland_FieldSelf_single.png",f)
    save("SingleIsland_out/SingleIsland_FieldSelf_single_poincare.png",h)
end




if save_reference_solution
    using JLD2
    save("SingleIsland_out/SingleIsland_FieldSelf_single.jld2","soln",soln)
end
=#

# q = Figure();
# ax_q = Axis(q[1,1]);
# s35 = lines!(ax_q,Dom1.gridx[:,1],soln35.u[2][:,floor(Int,ny/2)+1])
# s3 = lines!(ax_q,Dom1.gridx[:,1],soln3.u[2][:,floor(Int,ny/2)+1])
# s325 = lines!(ax_q,Dom1.gridx[:,1],soln325.u[2][:,floor(Int,ny/2)+1])
# s1 = scatter!(ax_q,Dom1.gridx[:,1],soln1.u[2][:,floor(Int,ny/2)+1])
# s0 = scatter!(ax_q,Dom1.gridx[:,1],soln0.u[2][:,floor(Int,ny/2)+1])
# s3n201 = lines!(ax_q,Dom.gridx[:,1],soln3n201.u[2][:,floor(Int,ny/2)+1])
# s35n201 = lines!(ax_q,Dom.gridx[:,1],soln35n201.u[2][:,floor(Int,ny/2)+1])

# Legend(q[1,2],[s35,s3,s325,s0,s3n201,s35n201],["3.5","3","3.25","0","3, n201","3.5, n201"])
# Legend(q[1,2],[s35,s3,s325,s0],["3.5","3","3.25","0"])

# surface!(ax_q,Dom.gridx,Dom.gridy,soln.u[2] .- soln2.u[2])

# f = Figure();
# ax_f = Axis(f[1,1]);
# ax2_f = Axis(f[1,2]);
# # surface!(ax_f,Dom.gridx,Dom.gridy,soln.u[2])
# lines!(ax_f,Dom.gridx[:,1],soln.u[2][:,floor(Int,ny/2)+1])
# lines!(ax2_f,Dom.gridx[:,1],soln2.u[2][:,floor(Int,ny/2)+1])
