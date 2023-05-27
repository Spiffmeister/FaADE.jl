
push!(LOAD_PATH,"../SBP_operators")
using SBP_operators


𝒟x = [0.0,1.0]
𝒟y = [-π,π]
nx = 21
ny = 21
Dom = Grid2D(𝒟x,𝒟y,nx,ny)

H = SBP_operators.Helpers.innerH(Dom,2)


u = rand(Dom.nx,Dom.ny);
v = rand(Dom.nx,Dom.ny);


H(u,v)

@code_warntype H(u,v)





χₘₙ = 2.1e-3 + 5.0e-3
params = (ϵₘₙ = [χₘₙ/2.0,χₘₙ/3.0 ],m=[2.0,3.0],n=[1.0,2.0])
function χ_h!(χ,x::Array{Float64},p,t)
    # Hamiltons equations for the field-line Hamiltonian
    # H = ψ²/2 - ∑ₘₙ ϵₘₙ(cos(mθ - nζ))
    χ[2] = x[1] #p_1            qdot        θ
    χ[1] = -sum(p.ϵₘₙ .*(sin.(p.m*x[2] - p.n*t) .* p.m)) #q_1        pdot        ψ
end

dH(X,x,p,t) = χ_h!(X,x,params,t)
PGrid = SBP_operators.construct_grid(dH,Dom,[-2π,2π])
Pfn = SBP_operators.generate_parallel_penalty(PGrid,Dom,2)






@code_warntype Pfn(u,v,0.1)




#=

=#
BoundaryLeft = Boundary(Dirichlet,(y,t) -> 0.0,Left,1)
BoundaryRight = Boundary(Dirichlet,(y,t) -> 1.0,Right,1)
BoundaryUpDown = PeriodicBoundary(2)
u₀(x,y) = x
order = 2

P = VariableCoefficientPDE2D(u₀,(x,y)->1e-8,(x,y)->1e-8,order,BoundaryLeft,BoundaryRight,BoundaryUpDown)

typeof(P.BoundaryConditions)
