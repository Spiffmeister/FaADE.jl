

abstract type newPDEProblem{dtype,DIMS} end





struct newProblem1D{TT      <: Real,
    DIMS,
    DCT     <: DiffusionCoefficient,
    ST      <: SourceTerm,
    DO      <: DerivativeOrder,
    SATB    <: SATBoundaries,
    PART    <: Union{ParallelData,Nothing}
        } <: newPDEProblem{TT,DIMS}
    InitialCondition    :: Function
    K                   :: DCT
    source              :: ST
    order               :: DO
    BoundaryConditions  :: SATB
    Parallel            :: PART

    function newProblem1D(order::Integer,u₀,K,G::GridType{TT,DIMS},BCs,S,Par) where {TT,DIMS}

        
        
        DiffCo = DiffusionCoefficient(K)
        source = SourceTerm{typeof(S)}(S)
        DiffOrd = DerivativeOrder{order}()

        # new{DIMS,TT}(u₀)
        new{TT,DIMS,typeof(DiffCo),typeof(source),typeof(DiffOrd),typeof(BCs),typeof(Par)}(u₀,DiffCo,source,DiffOrd,BCs,Par)
    end
end
newProblem1D(order,u₀,K,G,BCs) = newProblem1D(order,u₀,K,G,BCs,nothing,nothing)



struct newProblem2D{TT      <: Real,
        DIM,
        DCT     <: DiffusionCoefficient,
        ST      <: SourceTerm,
        DO      <: DerivativeOrder,
        SATB    <: SATBoundaries,
        PART    <: Union{ParallelData,Nothing}
            } <: newPDEProblem{TT,DIM}
    InitialCondition    :: Function
    Kx                  :: DCT
    Ky                  :: DCT
    source              :: ST
    order               :: DO
    BoundaryConditions  :: SATB
    Parallel            :: PART
    
    function newProblem2D(order::Integer,u₀,Kx,Ky,G::GridType{TT,DIM},BCs,S,Par) where {TT,DIM}

        source = SourceTerm{typeof(S)}(S)
        DO = DerivativeOrder{order}()
        DCx = DiffusionCoefficient(Kx)
        DCy = DiffusionCoefficient(Ky)


        # new(u₀,Kx,Ky,order,BoundaryConditions(Bounds))
        new{TT,2,typeof(DCx),typeof(source),typeof(DO),typeof(BCs),typeof(Par)}(u₀,DCx,DCy,source,DO,BCs)
    end
end




Base.ndims(P::newPDEProblem{TT,DIM}) where {TT,DIM} = DIM

Base.show(io::IO, P::newPDEProblem{TT,DIMS}) where {TT,DIMS} = print(io, DIMS," dimensional PDE Problem")



