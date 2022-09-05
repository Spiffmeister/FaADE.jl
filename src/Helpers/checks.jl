


@inline function check_order(order::Int)
    order ∈ [2,4,6] ? nothing : error("Order must be 2, 4 or 6")
end

@inline function check_boundary(side::NodeType)
    side ∈ [Left,Right] ? nothing : error("Boundary types must be Left or Right, cannot be Internal")
end

