
#=
struct FieldLineIntercept{F<:Union{Function,Nothing}}
    Intercept :: F
end


struct ParallelData{TT<:Real,
        DIM}
    PGrid       :: ParallelGrid
    κ           :: TT
    Intercept   :: FieldLineIntercept
    gridx       :: LinRange
    gridy       :: LinRange
    Δx          :: TT
    Δy          :: TT

    function ParallelData(PGrid::ParallelGrid,Grid::Grid2D{TT};κ=TT(1),intercept=nothing) where {TT}

        intercept_fieldlines = FieldLineIntercept(intercept)
        # intercept_fieldlines = intercept

        # K = DiffusionCoefficient(κ)

        gridx = LinRange(Grid.gridx[1],Grid.gridx[end],Grid.nx)
        gridy = LinRange(Grid.gridy[1],Grid.gridy[end],Grid.ny)

        new{TT,2}(PGrid,κ,intercept_fieldlines,gridx,gridy,Grid.Δx,Grid.Δy)
    end
end
=#

"""
    generate_parallel_penalty
Generates a 2D or 3D parallel penalty function given a parallel grid mapping.

Operator in the 1D case is found [here](@ref https://arxiv.org/abs/2303.15447). In the 2D case details can be found [here](@ref https://arxiv.org/abs/2306.00423).
"""
function generate_parallel_penalty end
function generate_parallel_penalty(planes::ParallelGrid{T,1,AT},order::Int;κ::T=1.0,τ::T=-1.0,interpmode::Symbol=:cust,interpfn::Union{Nothing,Function}=nothing) where {T,AT}
    # interp= choose_interpmode(interpmode=interpmode)
    if typeof(interpfn) == Nothing
        error("No interpolation function specified.")
    end

    H = build_H(order,grid.n)
    H = H.^-1.0
    
    α = -2.0
    τ = α/2.0

    let H=H, grid=grid, τ=τ, κ=κ, planes=planes

        ParPen(u,u₀,Δt) = ParallelPenalty1D!(interpfn,u,u₀,planes,Δt,grid,τ,κ,H)
    end
end
function generate_parallel_penalty(pargrid::ParallelGrid{T,2,AT},grid::Grid2D,order::Int;κ::T=1.0,τ::T=-1.0,perp::T=1.0,interpfn::Union{Nothing,Function}=nothing) where {T,AT}

    # interp = choose_interpmode(interpmode=interpmode)
    # if typeof(interpfn) == Nothing
    #     interp = LinearInterpolation
    # end

    H_x = build_H(order,grid.ny)
    H_y = build_H(order,grid.nx)

    H_x = 1.0 ./(H_x*grid.Δx)
    H_y = 1.0 ./(H_y*grid.Δy)

    # α = -2.0
    # τ = α/2.0 #TODO MIGHT BE A SQUARED HERE

    # itp(u) = LinearInterpolation((grid.gridx,grid.gridy),u)

    let H_x=H_x, H_y=H_y,
            τ=τ, κ=κ,
            nx=grid.nx, ny=grid.ny,
            D=grid,
            PGrid = pargrid,
            perp = perp
        ParPen(u,u₀,Δt) = ParallelPenalty2D!(u,u₀,Δt,PGrid,D,nx,ny,τ,κ,H_x,H_y,perp)
        # ParPen(u,u₀,Δt) = ParallelPenalty2D!(interp,u,u₀,Δt,PGrid,D,nx,ny,τ,κ,H_x,H_y)
        return ParPen
    end

end


"""
    ParallelPenalty1D
The default generated parallel penalty function for 2D (1D+parallel) problems
"""
function ParallelPenalty1D!(interp::Function,u::AbstractArray{TT},u₀::AbstractArray{TT},planes::ParallelGrid{TT,1},Δt::TT,grid::Grid1D,τ::TT,κ::TT,H::AbstractArray{TT}) where TT
    I = interp(grid.grid,u₀)
    for i = 1:grid.n
        u[i] = 1.0/(1.0 - κ*τ/2.0 * Δt * H[i]) * 
            (u[i] - Δt*κ*τ/4.0 * H[i] * (I(planes.FowardPlane[i]) + I(planes.BackwardPlane[i])))
    end
end
"""
    ParallelPenalty2D!(u::AbstractArray{T},u₀::AbstractArray{T},Δt::T,pargrid::ParallelGrid{2},grid::Grid2D,τ::T,κ::T,H_x::AbstractArray{T},H_y::AbstractArray{T}) where T
Default generated parallel penalty function for 3D (2D+parallel) problems
"""
function ParallelPenalty2D!(u::AbstractArray{TT},u₀::AbstractArray{TT},Δt::TT,
        PGrid::ParallelGrid{TT,2},D::Grid2D,
        nx::Integer,ny::Integer,τ::TT,κ::TT,H_x::AbstractArray{TT},H_y::AbstractArray{TT},perp::TT) where TT
    local wf :: TT
    local wb :: TT
    I = linear_interpolation((D.gridx,D.gridy),u)
    # I = scale(interpolate(u,BSpline(Linear())),(D.gridx,D.gridy))
    # I = scale(interpolate(u, BSpline(Quadratic(Line(OnGrid()))),),
        # range(D.gridx[1],D.gridx[end],length=D.nx),range(D.gridy[1],D.gridy[end],length=D.ny))
    # I = LinearInterpolation((D.gridx,D.gridy),u₀)

    # τ = -(1.e-14 + 1.0/κ)
    # τ = -(perp/κ)
    # τ = -sqrt(D.Δx * D.Δy/( (D.gridx[end]-D.gridx[1]) * (D.gridy[end] - D.gridy[1]) ))
    τ = -sqrt((D.gridx[end]-D.gridx[1])*(D.gridy[end]-D.gridy[1])/(D.Δx*D.Δy))  
    # -1.0/sqrt(D.Δx * D.Δy/( (D.gridx[end]-D.gridx[1]) * (D.gridy[end] - D.gridy[1]) ))
    # τ = -T(1)
    # τ = -D.Δx

    for j = 1:ny
        for i = 1:nx
            # H = H_y[i]*H_x[j]
            # H = 1.0
            wf = I(PGrid.Fplane.x[i,j],PGrid.Fplane.y[i,j])
            wb = I(PGrid.Bplane.x[i,j],PGrid.Bplane.y[i,j])
            u[i,j] = 1.0/(1.0 - κ* τ * Δt) * 
                ( u[i,j] -  κ*Δt*τ/2.0 * (wf + wb) )
            
            # u[i,j] = 1.0/(1.0 - κ* τ/2.0 * Δt * H) * 
            #     ( u[i,j] -  κ*Δt*τ/4.0 * H * 
            #         (I(PGrid.Fplane.x[i,j],PGrid.Fplane.y[i,j]) + I(PGrid.Bplane.x[i,j],PGrid.Bplane.y[i,j])) )
    
            
            # u[i,j] = u[i,j] + Δt * κ * τ/2.0 * H * (u[i,j] - 0.5*(I(PGrid.Fplane.x[i,j],PGrid.Fplane.y[i,j]) + I(PGrid.Bplane.x[i,j],PGrid.Bplane.y[i,j])))
            



        end
    end
end



"""
    For multiblock
"""
# function ParallelData()
# end
function applyParallelPenalty!(u::AbstractArray{TT},u₀::AbstractArray{TT},Δt::TT,θ::TT,
    P::ParallelData{TT,2}) where {TT}

    
    #= OLD METHOD =#
    #=
    local wf :: TT
    local wb :: TT

    I = scale( interpolate(u,BSpline(Cubic())),(P.gridx,P.gridy) )
    # I = scale( interpolate(u,BSpline(Quadratic())),(P.gridx,P.gridy) )
    # I = scale( interpolate(u,BSpline(Linear())),(P.gridx,P.gridy) )

    # τ = -sqrt((P.gridx[end]-P.gridx[1])*(P.gridy[end]-P.gridy[1])/(P.Δx*P.Δy))
    for j = 1:P.gridy.len
        for i = 1:P.gridx.len
            P.w_f[i,j] = I(P.PGrid.Fplane.x[i,j],P.PGrid.Fplane.y[i,j])
            P.w_b[i,j] = I(P.PGrid.Bplane.x[i,j],P.PGrid.Bplane.y[i,j])
        end
    end

    ττ = P.τ*norm(u - (P.w_f + P.w_b)/2,Inf)/norm((P.w_f + P.w_b)/2,Inf)^3

    @. u = 1.0/(1.0 - P.κ * ττ * Δt) * 
        ( u - P.κ*Δt*ττ/2.0 * (P.w_f + P.w_b) )
    =#


    # for j = 1:P.gridy.len
    #     for i = 1:P.gridx.len

    #         # wf = I(P.PGrid.Fplane.x[i,j],P.PGrid.Fplane.y[i,j])
    #         # wb = I(P.PGrid.Bplane.x[i,j],P.PGrid.Bplane.y[i,j])
    #         u[i,j] = 1.0/(1.0 - P.κ * ττ * Δt) * 
    #             ( u[i,j] -  P.κ*Δt*ττ/2.0 * (wf + wb) )
            
    #     end
    # end
    

    # w_f = getarray(P.w_f)
    # w_b = getarray(P.w_b)
    # println("parallel")
    #= NEW METHOD =#
    κ = P.κ
    w_b = P.w_b
    w_f = P.w_f
    H = P.H

    I = scale( interpolate(u,BSpline(Cubic())),(P.gridx,P.gridy) )
    # I = scale( interpolate(u,BSpline(Linear())),(P.gridx,P.gridy) )


    # w_f ← P_f u + P_b u
    _compute_w!(I,w_f,P.PGrid.Fplane,P.PGrid.Bplane,P.gridx.len,P.gridy.len)

    # Tune the parallel penalty
    τ = P.τ*(maximum(u - w_f)/maximum(w_f))^3
    # τ = 1e7*(maximum(abs.(u - w_f))/maximum(abs.(w_f)))^3
    # τ = P.τ
    P.τ_i[1] = τ
    # τ = P.τ


    # u^{n+1} = (1+θ)Δt κ τ
    # @. u = 1/(1 + θ * Δt * κ * τ * 1/H) * (u + θ*Δt*κ*τ/H * w_f) + (1-θ)*Δt * τ * κ * 1/H * (u₀ - w_b)
    # @show u[:,51]
    _compute_u!(u,w_f,w_b,κ,τ,θ,Δt,H,P.gridx.len,P.gridy.len)
    # @show u[:,51]

    # Replace the old parallel diffusion with the new one
    # TODO : FIX THIS FOR ADAPTIVE
    @. w_b = u - w_f

    # u = w_f
    u
end

function _compute_u!(u::AT,w_f::AT,w_b::AT,κ::TT,τ::TT,θ::TT,Δt::TT,H::CompositeH,nx::Int,ny::Int) where {TT,AT}
    # println(norm(u))
    for j in 1:ny
        for i in 1:nx
            u[i,j] = 1/(1 + θ * Δt * κ * τ / H[i,j]) * (
                (u[i,j] + θ * Δt * κ * τ / H[i,j] * w_f[i,j]) + 
                (1-θ)*Δt * κ * τ / H[i,j] * w_b[i,j] 
                )
        end
    end
    # println(norm(u))
    u
end

"""
    _compute_w!
Computes ``P_parallel u`` and stores it in `dest`.
"""
function _compute_w!(itp,dest::AT,Fplane::ParGrid,Bplane::ParGrid,nx::Int,ny::Int) where {AT}
    for j in 1:ny
        for i in 1:nx
            dest[i,j] = itp(Fplane.x[i,j],Fplane.y[i,j]) + itp(Bplane.x[i,j],Bplane.y[i,j])
            dest[i,j] = dest[i,j]/2
        end
    end
    dest
end

"""
    compute_parallel_operator
"""
function compute_parallel_operator(dest::AT,u::AT,P::ParallelData) where {AT}
    I = scale( interpolate(u,BSpline(Cubic())),(P.gridx,P.gridy) )
    _compute_w!(I,dest,P.PGrid.Fplane,P.PGrid.Bplane,P.gridx.len,P.gridy.len)
    @. dest = u - dest
end

# function P_parallel(dest,u,w,κ,τ)
#     @. dest = -κ * τ * P.H (u - w/2)
# end


