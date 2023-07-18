

abstract type newPDEProblem{DIMS,TT} end





struct newProblem1D{DIMS,
    TT      <: Real,
    GT      <: GridType,
    DCT     <: DiffusionCoefficient,
    ST      <: SourceTerm,
    DO      <: DerivativeOrder,
    SATB    <: SATBoundaries,
    PART    <: Union{ParallelData,Nothing}
        } <: newPDEProblem{DIMS,TT}
    InitialCondition    :: Function
    K                   :: DCT
    source              :: ST
    order               :: DO
    Grid                :: GT
    BoundaryConditions  :: SATB
    Parallel            :: PART

    function newProblem1D(order::Integer,u₀,K,G::GridType{TT,DIMS},BCs,S,Par) where {TT,DIMS}

        
        
        DiffCo = DiffusionCoefficient(K)
        source = SourceTerm{typeof(S)}(S)
        DiffOrd = DerivativeOrder{order}()

        # new{DIMS,TT}(u₀)
        new{DIMS,TT,typeof(G),typeof(DiffCo),typeof(source),typeof(DiffOrd),typeof(BCs),typeof(Par)}(u₀,DiffCo,source,DiffOrd,G,BCs,Par)
    end
end
newProblem1D(order,u₀,K,G,BCs) = newProblem1D(order,u₀,K,G,BCs,nothing,nothing)



struct newProblem2D{DIMS,
        TT      <: Real,
        GT      <: GridType,
        DCT     <: DiffusionCoefficient,
        ST      <: SourceTerm,
        DO      <: DerivativeOrder,
        SATB    <: SATBoundaries,
        PART    <: Union{ParallelData,Nothing}
            } <: newPDEProblem{DIMS,TT}
    InitialCondition    :: Function
    Kx                  :: DCT
    Ky                  :: DCT
    source              :: ST
    order               :: DO
    Grid                :: GT
    BoundaryConditions  :: SATB
    Parallel            :: PART
    
    function newProblem2D(order::Integer,u₀,Kx,Ky,G::GridType{TT,DIMS},BCs,S,Par) where {TT,KT<:Real,DIMS}

        source = SourceTerm{typeof(S)}(S)
        DO = DerivativeOrder{order}()
        DCx = DiffusionCoefficient(Kx)
        DCy = DiffusionCoefficient(Ky)


        # new(u₀,Kx,Ky,order,BoundaryConditions(Bounds))
        new{DIMS,TT,typeof(G),typeof(DCx),typeof(source),typeof(DO),typeof(BCs),typeof(Par)}(u₀,DCx,DCy,source,DO,G,BCs)
    end
end




Base.ndims(P::PDEProblem{D}) where D = D

Base.show(io::IO, P::newPDEProblem{DIMS,TT}) where {DIMS,TT} = print(io, DIMS," dimensional PDE Problem")



