

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

push!(LOAD_PATH,"."); using SBP_operators





𝒟 = [0.0,1.0]
n = 41
Dom = Grid1D(𝒟, n)


g₀(t) = 0.0
g₁(t) = 1.0
order = 2
K = ones(Float64,n)
Δt = 0.1Dom.Δx

# Define initial condition
u₀(x) = exp.(-(x.-0.5).^2 ./ 0.02)

# Define some boundary conditions
BoundaryDirichletLeft = Boundary(Dirichlet,g₀,Left,1)
BoundaryDirichletRight = Boundary(Dirichlet,g₁,Right,1)

# Build PDE problem
PD = VariableCoefficientPDE1D(u₀,K,order,BoundaryDirichletLeft,BoundaryDirichletRight)

# Testing internal data storage construction
BStor = SBP_operators.Helpers.BoundaryData1D{Float64}(PD.BoundaryConditions,order)
DStor = SBP_operators.Helpers.DataBlock{Float64}(PD.BoundaryConditions,Dom,Δt,2,PD.K)
CGStor = SBP_operators.Helpers.ConjGradBlock{Float64}(n)



# Testing internal boundary function construction
SATDL = SBP_operators.SATs.SAT_Dirichlet(PD.BoundaryConditions[1].RHS,Dom.Δx,PD.BoundaryConditions[1].side,PD.BoundaryConditions[1].axis,order)
SATDR = SBP_operators.SATs.SAT_Dirichlet(PD.BoundaryConditions[2].RHS,Dom.Δx,PD.BoundaryConditions[1].side,PD.BoundaryConditions[1].axis,order)
SATD1, SATFn1 = SBP_operators.SATs.SAT(PD.BoundaryConditions[1],Dom,order,:cgie)
SATD2, SATFn2 = SBP_operators.SATs.SAT(PD.BoundaryConditions[2],Dom,order,:cgie)



#= Neumann Boundaries =#
BoundaryNeumannLeft = Boundary(Neumann,g₀,Left,1)
BoundaryNeumannRight = Boundary(Neumann,g₁,Right,1)

PN = VariableCoefficientPDE1D(u₀,K,order,BoundaryNeumannLeft,BoundaryDirichletRight)

SATNL = SBP_operators.SATs.SAT_Neumann(PN.BoundaryConditions[1].RHS,Dom.Δx,PN.BoundaryConditions[1].side,PN.BoundaryConditions[1].axis,order)
SATNR = SBP_operators.SATs.SAT_Neumann(PN.BoundaryConditions[2].RHS,Dom.Δx,PN.BoundaryConditions[2].side,PN.BoundaryConditions[2].axis,order)
SATN1, SATNf1 = SBP_operators.SATs.SAT(PN.BoundaryConditions[1],Dom,order,:cgie)
SATN2, SATNf2 = SBP_operators.SATs.SAT(PN.BoundaryConditions[2],Dom,order,:cgie)







#####################


# Solution storage
IC = u₀(Dom.grid)
solution = SBP_operators.solvers.solution{Float64}(Dom,0.0,Δt,P)



function GenRHS(Dom,P)
        function RHS!(cache::AbstractArray,u::AbstractArray,k::AbstractArray) let n = Dom.n, Δx = Dom.Δx, o = P.order
                Dₓₓ!(cache,u,k,n,Δx,o)
                SATFn1(cache,u,k,SBP_operators.Helpers.SolutionMode)
                SATFn2(cache,u,k,SBP_operators.Helpers.SolutionMode)
            return nothing
        end
    end
end
function RHS!(cache::AbstractArray,u::AbstractArray,k::AbstractArray) let n = Dom.n, Δx = Dom.Δx, o = P.order
        Dₓₓ!(cache,u,k,n,Δx,o)
        SATFn1(cache,u,k,SBP_operators.Helpers.SolutionMode)
        SATFn2(cache,u,k,SBP_operators.Helpers.SolutionMode)
    return nothing
    end
end
rhs = GenRHS(Dom,P)

using BenchmarkTools

# @benchmark SBP_operators.solvers.A!(CGStor.rₖ,DStor.uₙ₊₁,RHS!,DStor.Δt,DStor.K);

@code_warntype RHS!(CGStor.rₖ,DStor.uₙ₊₁,DStor.K)
@code_warntype rhs!(CGStor.rₖ,DStor.uₙ₊₁,DStor.K)




CGTerm(cache::Array,u::Array,c::Array,::SBP_operators.Helpers.SATMode{:SolutionMode}) = 
            SAT_Dirichlet_implicit!(cache,SATDL.side,u,c,SATDL.α,SATDL.τ,SATDL.EDₓᵀ,P.order,eachcol)


CGTerm(cache::Array,u::Array,c::Array,::SBP_operators.Helpers.SATMode{:SolutionMode}) = 
            SAT_Dirichlet_implicit!(cache,SATDR.side,u,c,SATDR.α,SATDR.τ,SATDR.EDₓᵀ,P.order,eachcol)

@code_warntype CGTerm(DStor.uₙ₊₁,DStor.u,DStor.K,SBP_operators.Helpers.SolutionMode)






@time SBP_operators.solvers.A!(DStor.uₙ₊₁,DStor.u,RHS!,DStor.Δt,DStor.K)
@code_warntype SBP_operators.solvers.A!(DStor.uₙ₊₁,DStor.u,RHS!,DStor.Δt,DStor.K)


D = SBP_operators.Derivatives.generate_Derivative(41,0.25,2)

@code_warntype D(DStor.uₙ₊₁,DStor.u,DStor.K)

