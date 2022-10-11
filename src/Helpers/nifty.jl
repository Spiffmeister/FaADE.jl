#=
    NIFTY FUNCTIONS SUCH AS INPUT CHECKS, NODE ADJUSTMENTS ETC...
=#


#========
    CHECKING THINGS
========#
"""
    check_order
Ensure the users requested order makes sense
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
Returns the required number of nodes needed to compute `SecondDerivative` on the boundary
"""
@inline BoundaryNodeInput(order::Int) = order == 2 ? 1 : 2order
@inline BoundaryNodeInput(order::Int...) = BoundaryNodeInput.(order)

"""
    BoundaryNodeOutput
Returns the size of the array output by `SecondDerivative` on the boundary
"""
@inline BoundaryNodeOutput(order::Int) = order == 2 ? 1 : order+halforder(order)
@inline BoundaryNodeOutput(order::Int...) = BoundaryNodeOutput.(order)

"""
    SATNodeOutput
Returns the number of nodes needed for the `BoundaryData1D` and `BoundaryData2D` data structures.
"""
@inline SATNodeOutput(order::Int) = order == 2 ? 2 : order+halforder(order)
@inline SATNodeOutput(order::Int...) = SATNodeOutput.(order)


#========
    SELECT LOOPING DIRECTION
========#
function SelectLoopDirection(axis::Int)
    if axis == 1
        return eachcol
    elseif axis == 2
        return eachrow
    else
        error("axis must be 1 or 2")
    end
end


# @inline function Base.sum(u::AbstractArray{T},v::AbstractArray{T}) where T
#     tmp::T = 0
#     for (uᵢ,vᵢ) in zip(u,v)
#         tmp += uᵢ*vᵢ
#     end
# end

