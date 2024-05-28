

abstract type newPDEProblem{dtype,DIMS} end



struct BoundaryConditions{SATL<:SimultanousApproximationTerm,
        SATR<:Union{SimultanousApproximationTerm,Nothing},
        SATU<:Union{SimultanousApproximationTerm,Nothing},
        SATD<:Union{SimultanousApproximationTerm,Nothing}}
    BoundaryLeft    :: SATL
    BoundaryRight   :: SATR
    BoundaryUp      :: SATU
    BoundaryDown    :: SATD


    function BoundaryConditions(BCs...)
        if length(BCs) == 1
            new{Nothing,Nothing,Nothing,Nothing}(BCs,nothing,nothing,nothing)
        elseif length(BCs) == 2
            new{typeof(BCs[1]),typeof(BCs[2]),Nothing,Nothing}(BCs[1],BCs[2],nothing,nothing)
        elseif length(BCs) == 3
            new{typeof(BCs[1]),typeof(BCs[2]),typeof(BCs[3]),Nothing}(BCs[1],BCs[2],BCs[3],nothing)
        else
            new{typeof(BCs[1]),typeof(BCs[2]),typeof(BCs[3]),typeof(BCs[4])}(BCs[1],BCs[2],BCs[3],BCs[4])
        end
    end
end

"""
    SATBoundaries
"""
struct SATBoundaries{SATL<:SimultanousApproximationTerm,
        SATR<:Union{SimultanousApproximationTerm,Nothing},
        SATU<:Union{SimultanousApproximationTerm,Nothing},
        SATD<:Union{SimultanousApproximationTerm,Nothing}}
    BoundaryLeft    :: SATL
    BoundaryRight   :: SATR
    BoundaryUp      :: SATU
    BoundaryDown    :: SATD


    function SATBoundaries(BCs...)
        if length(BCs) == 1
            new{Nothing,Nothing,Nothing,Nothing}(BCs,nothing,nothing,nothing)
        elseif length(BCs) == 2
            new{typeof(BCs[1]),typeof(BCs[2]),Nothing,Nothing}(BCs[1],BCs[2],nothing,nothing)
        elseif length(BCs) == 3
            new{typeof(BCs[1]),typeof(BCs[2]),typeof(BCs[3]),Nothing}(BCs[1],BCs[2],BCs[3],nothing)
        else
            new{typeof(BCs[1]),typeof(BCs[2]),typeof(BCs[3]),typeof(BCs[4])}(BCs[1],BCs[2],BCs[3],BCs[4])
        end
    end
end













"""
    Problem1D
"""
struct Problem1D{TT      <: Real,
    DIMS,
    DCT,
    ST      <: SourceTerm,
    SATB    <: Union{Dict,Tuple},
    PART    # Parallel map, vector of parallel map, or nothing
        } <: newPDEProblem{TT,DIMS}
    InitialCondition    :: Function
    K                   :: DCT
    source              :: ST
    order               :: Int
    BoundaryConditions  :: SATB
    Parallel            :: PART

    function Problem1D(order::Integer,u₀,K,G::GridType{TT,DIMS},BCs,S,Par) where {TT,DIMS}

        source = SourceTerm{typeof(S)}(S)

        if typeof(G) <: GridMultiBlock
            typeof(BCs) <: Dict ? nothing : error("Boundary conditions must be a dictionary for multiblock problems")
            parse_boundaries(BCs,G)
        end

        # new{DIMS,TT}(u₀)
        new{TT,DIMS,typeof(K),typeof(source),typeof(BCs),typeof(Par)}(u₀,K,source,order,BCs,Par)
    end
end
Problem1D(order,u₀,K,G,BCs) = Problem1D(order,u₀,K,G,BCs,nothing,nothing)


"""
    Problem2D
"""
struct Problem2D{TT      <: Real,
        DIM,
        DCT,
        ST      <: SourceTerm,
        SATB    <: Union{Dict,Tuple},
        PART    # Parallel map, vector of parallel map, or nothing
            } <: newPDEProblem{TT,DIM}
    InitialCondition    :: Function
    Kx                  :: DCT
    Ky                  :: DCT
    source              :: ST
    order               :: Int
    BoundaryConditions  :: SATB
    Parallel            :: PART
    
    function Problem2D(order::Integer,u₀,Kx,Ky,G::GridType{TT,DIM},BCs,S,Par) where {TT,DIM}

        source = SourceTerm{typeof(S)}(S)

        if typeof(G) <: GridMultiBlock
            typeof(BCs) <: Dict ? nothing : error("Boundary conditions must be a dictionary for multiblock problems")
            parse_boundaries(BCs,G)
        end

        new{TT,2,typeof(Kx),typeof(source),typeof(BCs),typeof(Par)}(u₀,Kx,Ky,source,order,BCs,Par)
    end
end
Problem2D(order,u₀,Kx,Ky,G,BCs) = Problem2D(order,u₀,Kx,Ky,G,BCs,nothing,nothing)



Base.ndims(P::newPDEProblem{TT,DIM}) where {TT,DIM} = DIM

Base.show(io::IO, P::newPDEProblem{TT,DIMS}) where {TT,DIMS} = print(io, DIMS," dimensional PDE Problem")



