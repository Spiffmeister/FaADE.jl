
# abstract type DataBlockType{dtype<:AbstractFloat,N, atype<:AbstractArray{dtype}} end

abstract type LocalDataBlockType{dtype,DIM,atype} <: DataBlockType{dtype,DIM,atype} end


# abstract type Solution{T<:AbstractFloat} end


#========== WHOLE PROBLEM DATA ==========#
"""
    DataBlock
Passed around internally between functions. Only contains data required for current timestep.
"""
mutable struct DataBlock{T,DIM,AT,KT<:Union{AT,Vector{AT}}} <: LocalDataBlockType{T,DIM,AT}
    grid        :: GridType
    u           :: AT
    uₙ₊₁        :: AT
    K           :: KT
    boundary    :: BoundaryStorage{T,DIM,AT}
    t           :: T
    Δt          :: T
    Δu          :: T
    function DataBlock{T}(
            PDE::PDEProblem,
            grid::GridType{T,D},
            Δt::T) where {T,D}
    
        # If grid is 1D or 2D construct the right DataBlock
        if typeof(grid) <: Grid1D
            u   = zeros(T,grid.n)
            uₙ₊₁ = zeros(T,grid.n)
            BStor = BoundaryData1D{T}(PDE.BoundaryConditions,PDE.order)

            DiffCoeff = zeros(T,grid.n)
            DiffCoeff = PDE.K.(grid.grid)

            dim = 1
            
        elseif typeof(grid) <: Grid2D
            u   = zeros(T,(grid.nx,grid.ny))
            uₙ₊₁ = zeros(T,(grid.nx,grid.ny))
            BStor = BoundaryData2D{T}(PDE.BoundaryConditions,grid,PDE.order)

            DiffCoeff = [zeros(T,size(grid)),zeros(T,size(grid))]
            
            for j = 1:grid.ny
                for i = 1:grid.nx
                    DiffCoeff[1][i,j] = PDE.Kx(grid[i,j]...)
                    # DiffCoeff[1][i,j] = PDE.Kx(grid.gridx[i],grid.gridy[j])
                    DiffCoeff[2][i,j] = PDE.Ky(grid.gridx[i],grid.gridy[j])
                end
            end
            # for i in eachindex(grid)
            #     DiffCoeff[1][i] = PDE.Kx(grid[i]...)
            #     DiffCoeff[2][i] = PDE.Ky(grid[i]...)
            # end
            # setCoefficient!(PDE.Kx,DiffCoeff[1],grid)
            # setCoefficient!(PDE.Ky,DiffCoeff[2],grid)

            dim = 2

        end
        new{T,dim,typeof(u),typeof(DiffCoeff)}(grid,u,uₙ₊₁,DiffCoeff,BStor,0,Δt,0.0)
    end
end



#========== BOUNDARY DATA ==========#
"""
    BoundaryData1D
Data structure for storage of SATs in 1 dimensional problems
"""
mutable struct BoundaryData1D{T,AT} <: BoundaryStorage{T,1, AT}
    Type_Left   :: BoundaryConditionType
    Type_Right  :: BoundaryConditionType

    SAT_Left    :: AT
    SAT_Right   :: AT

    u_Left      :: AT
    u_Right     :: AT

    RHS_Left    :: AT
    RHS_Right   :: AT

    function BoundaryData1D{T}(BC::NamedTuple,order::Int) where {T}

        nnodes = SATNodeOutput(order)
        
        SAT_Left    = zeros(T,nnodes)
        SAT_Right   = zeros(T,nnodes)

        u_Left      = zeros(T,nnodes)
        u_Right     = zeros(T,nnodes)

        if length(BC) == 2
            BCL = BC.Left.type
            BCR = BC.Right.type
        elseif length(BC) == 1
            BCL = BCR = Periodic
        end

        new{T,typeof(u_Left)}(BCL,BCR,SAT_Left,SAT_Right,u_Left,u_Right,[0.0],[0.0])

    end
end

"""
    BoundaryData2D
Data structure for storage of SATs in 2 dimensional problems
"""
struct BoundaryData2D{T,AT} <: BoundaryStorage{T,2, AT}

    Type_Left    :: BoundaryConditionType
    Type_Right   :: BoundaryConditionType
    Type_Up      :: BoundaryConditionType
    Type_Down    :: BoundaryConditionType

    SAT_Left    :: AT #Same as u_ but for SAT
    SAT_Right   :: AT #Same as u_ but for SAT
    SAT_Up      :: AT #Same as u_ but for SAT
    SAT_Down    :: AT #Same as u_ but for SAT

    u_Left      :: AT #Solution along boundary with size determined by derivative order
    u_Right     :: AT #Solution along boundary with size determined by derivative order
    u_Up        :: AT #Solution along boundary with size determined by derivative order
    u_Down      :: AT #Solution along boundary with size determined by derivative order

    RHS_Left    :: AT
    RHS_Right   :: AT
    RHS_Up      :: AT
    RHS_Down    :: AT

    function BoundaryData2D{T}(BC::NamedTuple,grid::Grid2D,order::Int) where T

        nnodes = SATNodeOutput(order)
        ny = grid.ny
        nx = grid.nx

        SAT_Left =  zeros(T,(nnodes,ny)) #x
        SAT_Right = zeros(T,(nnodes,ny)) #x
        SAT_Up =    zeros(T,(nx,nnodes)) #y 
        SAT_Down =  zeros(T,(nx,nnodes)) #y

        u_Left =    zeros(T,(nnodes,ny)) #x
        u_Right =   zeros(T,(nnodes,ny)) #x
        u_Up =      zeros(T,(nx,nnodes)) #y 
        u_Down =    zeros(T,(nx,nnodes)) #y

        RHS_Left =    zeros(T,(1,ny)) #x
        RHS_Right =   zeros(T,(1,ny)) #x
        RHS_Up =      zeros(T,(nx,1)) #y 
        RHS_Down =    zeros(T,(nx,1)) #y

        new{T,typeof(u_Left)}(BC.Left.type,BC.Right.type,BC.Up.type,BC.Down.type,
            SAT_Left,SAT_Right,SAT_Up,SAT_Down,
            u_Left,u_Right,u_Up,u_Down,
            RHS_Left,RHS_Right,RHS_Up,RHS_Down)

    end
end




