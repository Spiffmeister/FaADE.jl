



struct Boundary
    type    :: BoundaryConditionType
    RHS     :: Function
    side    :: NodeType
    axis    :: Int
    function Boundary(type::BoundaryConditionType,RHS::Function,side::NodeType,axis::Int)

        type ∈ [Periodic,Interface] ? error("For periodic boundaries use PeriodicBoundary()") : nothing

        new(type,RHS,side,axis)
    end
end


struct PeriodicBoundary
    type    :: BoundaryConditionType
    axis    :: Int
    function PeriodicBoundary(axis::Int)
        if type != Periodic ? error("Periodic boundary must be Periodic.") : nothing

        new(Periodic,axis)
    end
end




struct problem
    PDE     :: Function
    SATs    :: Vector
end



