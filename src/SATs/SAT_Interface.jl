


"""
    SimultanousApproximationTerm
Abstract type for SAT data structures
"""
abstract type SimultanousApproximationTerm end



"""
    SAT(BoundCond::BoundaryConditionData,grid::GridType,order::Int,solver)
Creates a Dirichlet, Neumann or Periodic SAT function(s).

Inputs: 
- [`BoundaryConditionData`](@ref SPADE.Helpers.BoundaryConditionData)
- [`GridType`](@ref)
- Order of method
- Solver type (see [`solvers`](@ref))

Returns:
- BD: SAT Struct
- SATFn: SAT function(s)

If the solver is an implicit method, `SATFn` will have two methods, if explicit it will be a single method.
"""
function SAT(BoundCond::BoundaryConditionData,grid::GridType,order::Int,solver)
    # Is the problem 1 or 2D
    if typeof(grid) <: Grid1D
        Δ = grid.Δx
    elseif typeof(grid) <: Grid2D
        BoundCond.axis == 1 ? Δ=grid.Δx : Δ=grid.Δy
    end
    # Build Parametic things
    if BoundCond.type == Dirichlet
        BD = SAT_Dirichlet(BoundCond.RHS,Δ,BoundCond.side,BoundCond.axis,order)
    elseif BoundCond.type == Neumann
        BD = SAT_Neumann(BoundCond.RHS,Δ,BoundCond.side,BoundCond.axis,order)
    elseif BoundCond.type == Periodic
        BD = SAT_Periodic(Δ,BoundCond.axis,order)
    end
    # Build the SAT function
    SATFn = construct_SAT(BD,solver)
    return BD,SATFn
end

"""
    construct_SAT(Term::SimultanousApproximationTerm,solver)
Internal method for generating SAT functions.
"""
function construct_SAT(Term::SimultanousApproximationTerm,solver)
    if Term.type == Dirichlet
        SATFns = generate_Dirichlet(Term,solver)
    elseif Term.type == Neumann
        SATFns = generate_Neumann(Term,solver)
    elseif Term.type == Robin
        error("Robin not implemented")
    elseif Term.type == Periodic
        SATFns = generate_Periodic(Term,solver)
    else
        error("Available boundary types are Dirichlet, Neumann, Robin or Periodic")
    end
    return SATFns
end
















