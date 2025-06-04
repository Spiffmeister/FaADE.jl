using LinearAlgebra
using Printf
using Plots
using JLD2

using BenchmarkTools
using Profile
using PProf
using ProfileView
using Cthulhu

# cd("..")
using Interpolations
# push!(LOAD_PATH,"../FaADE")
using FaADE




###
𝒟x = [0.0,1.0]
𝒟y = [-π,π]
nx = 21
ny = 21
Dom = Grid2D(𝒟x,𝒟y,nx,ny)

kx(x,y) = 1.0
ky(x,y) = 1.0


Δt = 1.0 * min(Dom.Δx^2,Dom.Δy^2)
t_f = 10.0

u₀(x,y) = x


BoundaryLeft = Boundary(Dirichlet,(y,t) -> 0.0,Left,1)
BoundaryRight = Boundary(Dirichlet,(y,t) -> 1.0,Right,1)
BoundaryUpDown = PeriodicBoundary(2)


order = 2
method = :cgie


println(method)
println("(Δx,Δy)=(",Dom.Δx,",",Dom.Δy,")      ","Δt=",Δt,"        ","final time=",t_f)

P = VariableCoefficientPDE2D(u₀,kx,ky,order,BoundaryLeft,BoundaryRight,BoundaryUpDown)



χₘₙ = 2.1e-3 + 5.0e-3
params = (ϵₘₙ = [χₘₙ/2.0,χₘₙ/3.0 ],m=[2.0,3.0],n=[1.0,2.0])
function χ_h!(χ,x::Array{Float64},p,t)
    # Hamiltons equations for the field-line Hamiltonian
    # H = ψ²/2 - ∑ₘₙ ϵₘₙ(cos(mθ - nζ))
    χ[2] = x[1] #p_1            qdot        θ
    χ[1] = -sum(p.ϵₘₙ .*(sin.(p.m*x[2] - p.n*t) .* p.m)) #q_1        pdot        ψ
end

dH(X,x,p,t) = χ_h!(X,x,params,t)
PGrid = FaADE.construct_grid(dH,Dom,[-2π,2π])
Pfn = FaADE.generate_parallel_penalty(PGrid,Dom,2)


# using Profile
t_f = 1000Δt

# println("Benchmarking")
# @benchmark solve($P,$Dom,$Δt,$t_f,:cgie,penalty_func=$penalty_fn)
# using BenchmarkTools
# @time soln = solve(P,Dom,Δt,5.1Δt,:cgie,adaptive=true,penalty_func=penalty_fn)
# Profile.clear_malloc_data()

Pfn1 = FaADE.generate_parallel_penalty(PGrid,Dom,2)
P = VariableCoefficientPDE2D(u₀,(x,y)->1e-8,(x,y)->1e-8,order,BoundaryLeft,BoundaryRight,BoundaryUpDown)
# soln1 = solve(P,Dom,Δt,t_f,:cgie,adaptive=true,penalty_func=Pfn1)
PData = FaADE.ParallelData(PGrid,Dom)


# soln = solve(P,Dom,5.1Δt,t_f,:cgie,adaptive=true,Pgrid=PGrid)#,penalty_func=Pfn)
# Profile.clear_malloc_data()
# soln = solve(P,Dom,5.1Δt,t_f,:cgie,adaptive=true,Pgrid=PGrid)#,penalty_func=Pfn)



@profview soln1 = solve(P,Dom,Δt,2Δt,:cgie,adaptive=false,penalty_func=Pfn1)
@profview soln1 = solve(P,Dom,Δt,100.0,:cgie,adaptive=false,penalty_func=Pfn1)
@benchmark solve($P,$Dom,$Δt,$t_f,:cgie,adaptive=false,penalty_func=Pfn1)


@profview soln2 = solve(P,Dom,Δt,2Δt,:cgie,adaptive=false,Pgrid=PData)
@profview soln2 = solve(P,Dom,Δt,100.0,:cgie,adaptive=false,Pgrid=PData)
@benchmark solve($P,$Dom,$Δt,$t_f,:cgie,adaptive=false,Pgrid=PData)