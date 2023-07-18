#=
    NIFTY FUNCTIONS SUCH AS INPUT CHECKS, NODE ADJUSTMENTS ETC...
=#


#========
    CHECKING THINGS
========#
"""
    check_order
Ensure the users requested order makes sense.
"""
@inline function check_order(order::Int)
    order ∈ [2,4,6] ? nothing : error("Order must be 2, 4 or 6")
end
"""
    check_boundary
Some functions only work for `Left` or `Right` nodes.
"""
@inline function check_boundary(side::NodeType)
    side ∈ [Left,Right,Up,Down] ? nothing : error("Boundary types must be Left or Right, cannot be Internal")
end

"""
    check_boundary_args
Check the number of arguments supplied to a function.
"""
function check_boundary_args end
function check_boundary_args(f::Function)
    methods(f)[1].nargs - 1 == 2 ? nothing : error("Boundary conditions should be function with inputs f(x,t)")
end
function check_boundary_args(B::NamedTuple)
    for N in eachindex(B)
        if methods(B[N].RHS)[1].nargs - 1 == 2 
            error("Boundary conditions should be function with inputs f(x,t) at boundary ",N)
        end
    end
end


#========
    CHECKING THE NUMBER OF NODES REQUIRED GIVEN THE ORDER
========#
"""
    halforder
Returns ``order/2``
"""
@inline halforder(order::Int) = Int(order/2)
@inline halforder(order::Int...) = halforder.(order)

"""
    BoundaryNodeInput
Returns the required number of nodes needed to compute [`FaADE.Derivatives.SecondDerivative`](@ref) on the boundary.
``order == 2, return 1``, ``order > 2, return 2order``
"""
@inline BoundaryNodeInput(order::Int) = order == 2 ? 1 : 2order
@inline BoundaryNodeInput(order::Int...) = BoundaryNodeInput.(order)

"""
    BoundaryNodeOutput
Returns the size of the array output by [`FaADE.Derivatives.SecondDerivative`](@ref) on the boundary.
``order == 2, return 1``, ``order > 2, return order+order/2``
"""
@inline BoundaryNodeOutput(order::Int) = order == 2 ? 1 : order+halforder(order)
@inline BoundaryNodeOutput(order::Int...) = BoundaryNodeOutput.(order)

"""
    SATNodeOutput
Returns the number of nodes needed for the [`BoundaryData1D`](@ref) and [`BoundaryData2D`](@ref) data structures.
``order == 2, return 2``, ``order > 2, return order+order/2``
"""
@inline SATNodeOutput(order::Int) = order == 2 ? 2 : order+halforder(order)
@inline SATNodeOutput(order::Int...) = SATNodeOutput.(order)


#========
    SELECT LOOPING DIRECTION
========#
"""
    SelectLoopDirection
Used to tell the SATs over which axis to loop
"""
function SelectLoopDirection(axis::Int)
    if axis == 1
        return eachcol
    elseif axis == 2
        return eachrow
    else
        error("axis must be 1 or 2")
    end
end


"""
    GetAxis
Return the axis along which the node lies
"""
GetAxis(::NodeType{T,Ax}) where {T,Ax} = Int(Ax)

"""
    GetDim
Return the data structures associated dimension
"""
function GetDim end
GetDim(::BoundaryStorage{T,D,AT}) where {T,D,AT} = Int(D)
GetDim(::DataBlockType{T,D,AT}) where {T,D,AT} = Int(D)
