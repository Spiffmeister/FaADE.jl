
"""
    PDEProblem{dtype,DIMS}
Abstract type for `PDEProblemND` structs.

see [`Problem1D`](@ref) and [`Problem2D`](@ref)
"""
abstract type PDEProblem{dtype,DIMS} end





"""
    Problem1D{TT      <: Real,
    DIMS,
    DCT,
    ST      <: SourceTerm,
    SATB    <: Union{Dict,Tuple},
    PART} <: PDEProblem{TT,DIMS}

Struct for one dimensional problems
"""
struct Problem1D{TT      <: Real,
    DIMS,
    DCT,
    ST      <: SourceTerm,
    SATB    <: Union{Dict,Tuple},
    PART    # Parallel map, vector of parallel map, or nothing
        } <: PDEProblem{TT,DIMS}
    InitialCondition    :: Function
    K                   :: DCT
    source              :: ST
    order               :: Int
    BoundaryConditions  :: SATB
    Parallel            :: PART

    @doc """
        Problem1D(order::Integer,u₀,K,G::GridType{TT,DIMS},BCs,S,Par) where {TT,DIMS}
    
    Required arguments:
    - `order`: order of the solver to be used
    - `u₀`: Initial condition
    - `K`: Diffusion coefficient
    - `G`: one dimensional [`Grid1D`](@ref) or [`GridMultiBlock`](@ref)
    - `BCs`: Tuple of boundary conditions (see [`SATs`](@ref))
    """
    function Problem1D(order::Integer,u₀,K,G::GridType{TT,DIMS},BCs;source=nothing,parallel=nothing) where {TT,DIMS}

        S = SourceTerm{typeof(source)}(source)

        if typeof(G) <: GridMultiBlock
            typeof(BCs) <: Dict ? nothing : error("Boundary conditions must be a dictionary for multiblock problems")
            parse_boundaries(BCs,G)
        end

        # new{DIMS,TT}(u₀)
        new{TT,DIMS,typeof(K),typeof(S),typeof(BCs),typeof(parallel)}(u₀,K,S,order,BCs,parallel)
    end
end


"""
    Problem2D{TT      <: Real,
        DIM,
        DCT,
        ST      <: SourceTerm,
        SATB    <: Union{Dict,Tuple},
        PART} <: PDEProblem{TT,DIM}

Struct for two dimensional problems
"""
struct Problem2D{TT      <: Real,
        DIM,
        DCT,
        ST      <: SourceTerm,
        SATB    <: Union{Dict,Tuple},
        PART    # Parallel map, vector of parallel map, or nothing
            } <: PDEProblem{TT,DIM}
    InitialCondition    :: Function
    Kx                  :: DCT
    Ky                  :: DCT
    source              :: ST
    order               :: Int
    BoundaryConditions  :: SATB
    Parallel            :: PART
    
    @doc """
        Problem2D(order::Integer,u₀,Kx,Ky,G::GridType{TT,DIM},BCs,S,Par) where {TT,DIM}
    
    Required arguments:
    - `order`: order of the solver to be used
    - `u₀`: Initial condition
    - `Kx`: Diffusion coefficient in ``x``
    - `Ky`: Diffusion coefficient in ``y``
    - `G`: one dimensional [`Grid2D`](@ref) or [`GridMultiBlock`](@ref)
    - `BCs`: Tuple of boundary conditions (see [`SATs`](@ref))
    """
    function Problem2D(order::Integer,u₀,Kx,Ky,G::GridType{TT,DIM},BCs;source=nothing,parallel=nothing) where {TT,DIM}

        S = SourceTerm{typeof(source)}(source)

        if typeof(G) <: GridMultiBlock
            typeof(BCs) <: Dict ? nothing : error("Boundary conditions must be a dictionary for multiblock problems")
            parse_boundaries(BCs,G)
        end

        new{TT,2,typeof(Kx),typeof(S),typeof(BCs),typeof(parallel)}(u₀,Kx,Ky,S,order,BCs,parallel)
    end
end



Base.ndims(P::PDEProblem{TT,DIM}) where {TT,DIM} = DIM

Base.show(io::IO, P::PDEProblem{TT,DIMS}) where {TT,DIMS} = print(io, DIMS," dimensional PDE Problem")



