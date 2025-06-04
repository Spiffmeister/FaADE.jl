using LinearAlgebra
using Revise

using FaADE


order = 2

K = 1.0
n = 41

ωx = 1.0
ωt = 1.0
cx = 0.0

K = 1.0

θ = 0.5


# Solution

exact(x,t) = cos(2π*ωt*t) * sin(2π*x*ωx + cx)
u₀(x) = exact(x,0.0)
F(x,t) = -2π*ωt*sin(2π*ωt*t)*sin(2π*x*ωx + cx) + K * 4π^2 * ωx^2 * cos(2π*ωt*t)*sin(2π*x*ωx + cx)
# DIRICHLET
BxL(t) = cos(2π*ωt*t) * sin(cx) #Boundary condition x=0
BxR(t) = cos(2π*ωt*t) * sin(2π*ωx + cx) #Boundary condition x=Lx
# NEUMANN
# BxL(t) = 2π*ωx * K * cos(2π*ωt*t) * cos(cx) #Boundary condition x=0
# BxR(t) = 2π*ωx * K * cos(2π*ωt*t) * cos(2π*ωx + cx) #Boundary condition x=Lx




#====== New solver 1 volume ======#
Dom = Grid1D([0.0,1.0],n)
Dl = FaADE.SATs.SAT_Dirichlet(BxL,Dom.Δx,Left,  order)
Dr = FaADE.SATs.SAT_Dirichlet(BxR,Dom.Δx,Right, order)
# BD = (Dl,Dr)
BD = (Dl,Dr)

P = Problem1D(order,u₀,K,Dom,BD,F,nothing)
println("---Solving 1 volume---")
soln = solve(P,Dom,Δt,t,solver=:theta,θ=θ)


#====== New solver 2 volume ======#
println("2 volume")
# D1 = Grid1D([0.0,0.5],floor(Int,n/2)+1)
# D2 = Grid1D([0.5,1.0],floor(Int,n/2)+1)
D1 = Grid1D([0.0,0.5],61)
D2 = Grid1D([0.5,1.0],31)

Dom2V = GridMultiBlock(D1,D2)

Dl = FaADE.SATs.SAT_Dirichlet(BxL,D1.Δx,Left, order)
Dr = FaADE.SATs.SAT_Dirichlet(BxR,D2.Δx,Right,order)

BD2 = Dict(1 => (Dl,), 2 => (Dr,))

# BD = (Dl,Dr)
println("---Solving 2 volume---")
P2V = Problem1D(order,u₀,K,Dom2V,BD2,F,nothing)

soln2V = solve(P2V,Dom2V,Δt,t,solver=:theta,θ=θ)


#====== New solver 3 volume ======#
println("3 volume")
D1 = Grid1D([0.0,0.40],51)
D2 = Grid1D([0.40,0.60],41)
D3 = Grid1D([0.60,1.0],51)


Dom3V = GridMultiBlock(D1,D2,D3)

Dl = FaADE.SATs.SAT_Dirichlet(BxL,D1.Δx,Left,    order)
Dr = FaADE.SATs.SAT_Dirichlet(BxR,D3.Δx,Right,   order)
# BD = (Dl,Dr)

BD3 = Dict(1 => (Dl,), 3 => (Dr,))

P3V = Problem1D(order,u₀,K,Dom3V,BD3,F,nothing)

println("---Solving 3 volume---")
soln3V = solve(P3V,Dom3V,Δt,t)



e = [exact(Dom[i],t) for i in eachindex(Dom)]
u0 = [u₀(Dom[i]) for i in eachindex(Dom)]


