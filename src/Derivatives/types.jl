

"""
Abstract type for derivative operators
"""
abstract type DerivativeOperatorType{DIM} end



"""
    DiffusionOperator{TT<:Real,DO,COEFF} <: DerivativeOperatorType{1}

Data structure for diffusion operators.
"""
struct DiffusionOperator{TT<:Real,DO,COEFF} <: DerivativeOperatorType{1}
    order   :: Int64
    n       :: Int64
    Δx      :: TT
    periodic:: Bool

    @doc """
        DiffusionOperator(n::Int64,Δx::TT,order::Int,periodic::Bool,coeff::Symbol) where {TT<:Real}

    Required arguments:
    - `n`: Grid size being applied to
    - `\\Delta x`: Grid spacing
    - `order`: SBP operator order
    - `periodic`: Is this a periodic stencil
    - `coeff`: `:constant` or `:variable` coefficient operator
    """
    function DiffusionOperator(n::Int64,Δx::TT,order::Int,periodic::Bool,coeff::Symbol) where {TT<:Real}
        return new{TT,order,coeff}(order,n,Δx,periodic)
    end
end
"""
    DiffusionOperatorND{TT<:Real,DIM,DO,COEFF,AT<:AbstractArray{TT,DIM}} <: DerivativeOperatorType{DIM}

N-dimensional diffusion operator.
"""
struct DiffusionOperatorND{TT<:Real,DIM,DO,COEFF,AT<:AbstractArray{TT,DIM}} <: DerivativeOperatorType{DIM}
    DO      :: NTuple{DIM,DiffusionOperator{TT,DO,COEFF}}
    cache   :: AT
    @doc """
        DiffusionOperatorND(D::DiffusionOperator{TT,DO,COEFF}...) where {TT<:Real,DO,COEFF}
    
    Input is

    ```
        julia> Dxx = DiffusionOperator(...)
        julia> Dyy = DiffusionOperator(...)

        julia> D = DiffusionOperatorND(Dxx,Dyy)
    ```
    """
    function DiffusionOperatorND(D::DiffusionOperator{TT,DO,COEFF}...) where {TT<:Real,DO,COEFF}
        DIM = length(D)
        # cache = zeros(AT,(DO[1].n,DO[2].n))
        if COEFF == :Constant
            cache = zeros(TT,(1,1))
        elseif COEFF == :Variable
            cache = zeros(TT,(D[1].n,D[2].n))
        end
        
        new{TT,DIM,D[1].order,COEFF,typeof(cache)}(D,cache)
    end
end

# Base.show(io::IO,DO::DerivativeOperator{TT,DIM,O}) where {TT,DIM,O} = print("Order ",O," ",DIM," dimensional diffusion SBP operator.")
Base.show(io::IO,DO::DiffusionOperator{TT,O,COEFF}) where {TT,O,COEFF} = print("Order ",O," diffusion SBP operator.")



# _next(i,start,stop) = i+1 > stop ? rem(i+1,stop+1)+start : i+1
# _prev(i,stop) = i-1 ≤ start ? mod1(i-1,stop) : i-1
@inline _next(i,stop) = i+1 > stop ? rem(i+1,stop+1)+2 : i+1
@inline _prev(i,stop) = i-1 ≤ 1 ? mod1(i-1,stop-1) : i-1

# Base.show(io::IO,DO::DerivativeOrder{O}) where O = print(io, "order ",O," second derivative operator")
