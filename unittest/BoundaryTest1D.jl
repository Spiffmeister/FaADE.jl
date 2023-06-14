

#====
Testing the boundary operator construction

For Problem setup:
    1. Select PDE type (only handles variable coefficient diffusion problem at the moment)
        a. Order
        b. Domain
        c. Grid

For Dirichlet, Neumann, Robin the workflow is:
    1. Set 
====#

push!(LOAD_PATH,"."); using SPADE



𝒟 = [0.0,1.0]
n = 41
Dom = Grid1D(𝒟, n)

order = 2
K = ones(Float64,n)
Δt = 0.1Dom.Δx

# Define initial condition
u₀(x) = exp.(-(x.-0.5).^2 ./ 0.02)

# Define some boundary conditions
BoundaryDirichletLeft = Boundary(Dirichlet,t->0.0,Left,1)
BoundaryDirichletRight = Boundary(Dirichlet,t->1.0,Right,1)

# Build PDE problem
PD = VariableCoefficientPDE1D(u₀,K,order,BoundaryDirichletLeft,BoundaryDirichletRight)

# Testing internal data storage construction
BStor = SPADE.Helpers.BoundaryData1D{Float64}(PD.BoundaryConditions,order)
DStor = SPADE.Helpers.DataBlock{Float64}(PD.BoundaryConditions,Dom,Δt,2,PD.K)
CGStor = SPADE.Helpers.ConjGradBlock{Float64}(Dom,order)

x = rand(n); y = ones(n);

CGStor.innerprod(x,y) == (x[1]*0.5 + x[end]*0.5 + sum(x[2:end-1])) * Dom.Δx



# Testing internal boundary function construction
SATDL = SPADE.SATs.SAT_Dirichlet(PD.BoundaryConditions[1].RHS,Dom.Δx,PD.BoundaryConditions[1].side,PD.BoundaryConditions[1].axis,order)
SATDR = SPADE.SATs.SAT_Dirichlet(PD.BoundaryConditions[2].RHS,Dom.Δx,PD.BoundaryConditions[2].side,PD.BoundaryConditions[1].axis,order)
SATD1, SATFn1 = SPADE.SATs.SAT(PD.BoundaryConditions[1],Dom,order,:cgie)
SATD2, SATFn2 = SPADE.SATs.SAT(PD.BoundaryConditions[2],Dom,order,:cgie)



#= Neumann Boundaries =#
BoundaryNeumannLeft = Boundary(Neumann,g₀,Left,1)
BoundaryNeumannRight = Boundary(Neumann,g₁,Right,1)

PN = VariableCoefficientPDE1D(u₀,K,order,BoundaryNeumannLeft,BoundaryDirichletRight)

SATNL = SPADE.SATs.SAT_Neumann(PN.BoundaryConditions[1].RHS,Dom.Δx,PN.BoundaryConditions[1].side,PN.BoundaryConditions[1].axis,order)
SATNR = SPADE.SATs.SAT_Neumann(PN.BoundaryConditions[2].RHS,Dom.Δx,PN.BoundaryConditions[2].side,PN.BoundaryConditions[2].axis,order)
SATN1, SATNf1 = SPADE.SATs.SAT(PN.BoundaryConditions[1],Dom,order,:cgie)
SATN2, SATNf2 = SPADE.SATs.SAT(PN.BoundaryConditions[2],Dom,order,:cgie)




#= Periodic Boundaries =#
BoundaryPeriodicTest = PeriodicBoundary(1)

PP = VariableCoefficientPDE1D(u₀,K,order,BoundaryPeriodicTest)

SATP = SPADE.SATs.SAT_Periodic(Dom.Δx,PP.BoundaryConditions[1].axis,order)

SATPf1 = SPADE.SATs.SAT(PP.BoundaryConditions[1],Dom,order,:cgie)




#####################


# Solution storage
IC = u₀(Dom.grid)
solution = SPADE.solvers.solution{Float64}(Dom,0.0,Δt,P)



function GenRHS(Dom,P)
        function RHS!(cache::AbstractArray,u::AbstractArray,k::AbstractArray) let n = Dom.n, Δx = Dom.Δx, o = P.order
                D₂!(cache,u,k,n,Δx,o)
                SATFn1(cache,u,k,SPADE.Helpers.SolutionMode)
                SATFn2(cache,u,k,SPADE.Helpers.SolutionMode)
            return nothing
        end
    end
end
function RHS!(cache::AbstractArray,u::AbstractArray,k::AbstractArray) let n = Dom.n, Δx = Dom.Δx, o = P.order
        D₂!(cache,u,k,n,Δx,o)
        SATFn1(cache,u,k,SPADE.Helpers.SolutionMode)
        SATFn2(cache,u,k,SPADE.Helpers.SolutionMode)
    return nothing
    end
end
rhs = GenRHS(Dom,P)

using BenchmarkTools

# @benchmark SPADE.solvers.A!(CGStor.rₖ,DStor.uₙ₊₁,RHS!,DStor.Δt,DStor.K);

@code_warntype RHS!(CGStor.rₖ,DStor.uₙ₊₁,DStor.K)
@code_warntype rhs!(CGStor.rₖ,DStor.uₙ₊₁,DStor.K)




CGTerm(cache::Array,u::Array,c::Array,::SPADE.Helpers.SATMode{:SolutionMode}) = 
            SAT_Dirichlet_implicit!(cache,SATDL.side,u,c,SATDL.α,SATDL.τ,SATDL.ED₁ᵀ,P.order,eachcol)


CGTerm(cache::Array,u::Array,c::Array,::SPADE.Helpers.SATMode{:SolutionMode}) = 
            SAT_Dirichlet_implicit!(cache,SATDR.side,u,c,SATDR.α,SATDR.τ,SATDR.ED₁ᵀ,P.order,eachcol)

@code_warntype CGTerm(DStor.uₙ₊₁,DStor.u,DStor.K,SPADE.Helpers.SolutionMode)






@time SPADE.solvers.A!(DStor.uₙ₊₁,DStor.u,RHS!,DStor.Δt,DStor.K)
@code_warntype SPADE.solvers.A!(DStor.uₙ₊₁,DStor.u,RHS!,DStor.Δt,DStor.K)


D = SPADE.Derivatives.generate_Derivative(41,0.25,2)

@code_warntype D(DStor.uₙ₊₁,DStor.u,DStor.K)

