"""
    applyParallelPenalty!
The default generated parallel penalty function for 2D (1D+parallel) problems

Operator in the 1D case is found [here](https://arxiv.org/abs/2303.15447). In the 2D case details can be found [here](https://arxiv.org/abs/2306.00423).    
"""
function applyParallelPenalty! end

"""
1D apply parallel penalty
"""
function applyParallelPenalty!(u::VT,t::TT,Δt::TT,
        P::ParallelData{TT,1},grid::Grid1D{TT,MET}) where {TT,VT,MET}
    
    κ = P.κ
    w = P.w
    H = P.H

    I = LinearInterpolator(P.gridx,u)

    Fplane = P.PGrid.Fplane
    Bplane = P.PGrid.Bplane

    # w ← P_f u + P_b u
    for i in eachindex(P.gridx)
        x_f = Fplane.x[i]
        x_b = Bplane.x[i]
        if grid[1] ≤ x_f ≤ grid[grid.n]
            w[i] = I(x_f)
        else
            w[i] = TT(0)
        end
        if grid[1] ≤ x_b ≤ grid[grid.n]
            w[i] += I(x_b)
        else
            w[i] += TT(0)
        end
    end

    τ = P.τ/grid.Δx # Tune the parallel penalty

    P.w_b[1,1] = maximum(u - w)

    for i in 1:grid.n
        u[i] = TT(1)/(TT(1) - τ*κ/TT(2) * Δt * H.H[1][i]) *
            (u[i] - Δt * κ * τ / TT(4) * H.H[1][i] * w[i])
    end

    u
end
"""
2D Parallel penalty for single block problems
"""
function applyParallelPenalty!(u::AbstractArray{TT},t::TT,Δt::TT,
        P::ParallelData{TT,2,PGT,GT,BT,IT},grid::Grid2D{TT,MET}) where {TT,MET,PGT,GT,BT,IT}

    κ = P.κ
    w = P.w
    H = P.H
    J = grid.J

    I = P.Interpolant
    IC = P.Intercept

    # w ← P_f u + P_b u
    if isnothing(IC)
        _compute_w!(I,w,P.PGrid.Fplane,P.PGrid.Bplane,grid.nx,grid.ny)
    else
        _compute_w!(I,IC,w,P.PGrid.Fplane,P.PGrid.Bplane,t,grid.nx,grid.ny)
    end
    
    # Tune the parallel penalty
    τ = P.τ*0.1*(maximum(abs.(u - w))/ maximum(abs.(w)))^2.0
    isinf(τ) ? τ = TT(1) : nothing
    P.τ_i[1] = τ * κ * 1/(maximum(J) * P.Δx*P.Δy * maximum(H.H[1].Boundary) * maximum(H.H[2].Boundary))

    # P.w_b[1,1] = maximum(u - w)

    # u^{n+1} = (1+θ)Δt κ τ
    if MET == CartesianMetric
        _compute_u!(u,w,κ,τ,Δt,H,grid.nx,grid.ny)
    else
        _compute_u_curvilinear!(u,w,κ,τ,Δt,H,grid.nx,grid.ny,J)
    end

    u
end
"""
Applies the parallel penalty for a multiblock problem
"""
function applyParallelPenalty!(u::AbstractArray{TT},τ::TT,Δt::TT,
        P::Dict{Int64,ParallelData{TT,DIM}},grid::Grid2D{TT,MET},I) where {TT,MET,DIM}

    LocalP = P[I] :: ParallelData{TT,DIM}

    κ = LocalP.κ
    w = LocalP.w
    H = LocalP.H
    Jac = grid.J

    ny = grid.ny :: Integer
    nx = grid.nx :: Integer

    if length(Jac) > 1
        for j in Base.OneTo(ny)#1:grid.ny
            for i in Base.OneTo(nx)#1:grid.nx
                u[i,j] = 1/(1 + Δt * κ * τ / (Jac[i,j] * H[i,j])) * (
                    u[i,j] + Δt * κ * τ * w[i,j] / (Jac[i,j] * H[i,j])
                )
            end
        end
    else
        for j in Base.OneTo(ny)#1:grid.ny
            for i in Base.OneTo(nx)#1:grid.nx
                u[i,j] = 1/(1 + Δt * κ * τ / H[i,j]) * (
                    u[i,j] + Δt * κ * τ * w[i,j] / H[i,j]
                )
            end
        end
    end

    # τ = P[I].τ * 0.1 * (maximum(abs.(u - w))/ maximum(abs.(w)))^2.0

end

function _compute_u!(u::AT,w::AT,κ::TT,τ::TT,Δt::TT,H::CompositeH,nx::Int,ny::Int) where {TT,AT}
    # println(norm(u))
    # um = maximum(abs.(u .- w))
    for j in 1:ny
        for i in 1:nx
            ### CURRENTLY IN THE PAPER

            u[i,j] = 1/(1 + Δt * κ * τ / (H[i,j])) * (
                u[i,j] + 
                +Δt*κ*τ*w[i,j]/(H[i,j])
                    )

        end
    end
    u
end

function _compute_u_curvilinear!(u::AT,w::AT,κ::TT,τ::TT,Δt::TT,H::CompositeH,nx::Int,ny::Int,J::AT) where {TT,AT}
    for j in 1:ny
        for i in 1:nx

            u[i,j] = 1/(1 + Δt * κ * τ / (J[i,j]*H[i,j])) * (
                u[i,j] + 
                +Δt*κ*τ*w[i,j]/(J[i,j]*H[i,j])
                )
        end
    end
    u
end

"""
    _compute_w!
Computes ``P_{parallel} u`` and stores it in `dest`.
"""
function _compute_w! end
"""
Compute `w = (w + w_b)/2` when an interpolating function is provided.
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
If an `intercept(x,y,t)` function is provided then call it each iteration. To be safe ensure `intercept(u,x)` returns `u` if the intercept has no effect.
"""
function _compute_w!(itp,itc,dest::AT,Fplane::ParGrid,Bplane::ParGrid,t::TT,nx::Int,ny::Int) where {AT,TT}
    for j in 1:ny
        for i in 1:nx
            dest[i,j] = itp(Fplane.x[i,j],Fplane.y[i,j]) + itp(Bplane.x[i,j],Bplane.y[i,j])
            dest[i,j] = itc(dest[i,j],Fplane.x[i,j],Fplane.y[i,j],t) + itc(dest[i,j],Bplane.x[i,j],Bplane.y[i,j],t)
            dest[i,j] = dest[i,j]/2
        end
    end
    dest
end
function _compute_w!(itp,dest::AT,plane::ParGrid,nx::Int,ny::Int) where {AT}
    for j in 1:ny
        for i in 1:nx
            dest[i,j] = itp(plane.x[i,j],plane.y[i,j])
        end
    end
    dest
end
function _compute_w!(itp,dest::AT,plane::ParGrid) where {AT}
    for j in 1:plane.y.len
        for i in 1:plane.x.len
            dest[i,j] = itp(plane.x[i,j],plane.y[i,j])
        end
    end
    dest
end


"""
    compute_parallel_operator
"""
function compute_parallel_operator(dest::AT,u::AT,P::ParallelData) where {AT}
    I = BicubicInterpolator(P.gridx,P.gridy,u)
    _compute_w!(I,dest,P.PGrid.Fplane,P.PGrid.Bplane,length(P.gridx),length(P.gridy))
    @. dest = u - dest
end



"""
    computeglobalw!
"""
function computeglobalw! end

function computeglobalw!(PD::ParallelMultiBlock{TT,DIM,IT,CT},uglobal::VAT,t::TT,Δt::TT) where {TT,VAT,DIM,IT,CT}

    Interpolant = PD.Interpolant
    Intercept = PD.Intercept
    τglobal = PD.τ

    Threads.@threads for I in eachindex(Interpolant)
        # TODO: Place a function barrier here probably, there is a type instability
        #   caused by the selection of localP
        localP = PD.PData[I] :: ParallelData{TT,DIM}
        w = localP.w
        PGrid = localP.PGrid :: ParallelGrid
        computeglobalw!(Interpolant,Intercept,w,uglobal[I],PGrid,τglobal,t,I)
    end
    τglobal .= τglobal * PD.PData[1].τ
end
function computeglobalw!(interpolant::IT,intercept::CT,w::AT,u::AT,PGrid::ParallelGrid,τglobal::VT,t::TT,I::Int) where {AT,VT,TT,IT,CT}

    sgiF = PGrid.Fplane.subgrid
    sgiB = PGrid.Bplane.subgrid

    nnFx = PGrid.Fplane.x
    nnFy = PGrid.Fplane.y
    nnBx = PGrid.Bplane.x
    nnBy = PGrid.Bplane.y

    for J in eachindex(w)
        tmpf = interpolant[sgiF[J]](nnFx[J],nnFy[J])
        tmpb = interpolant[sgiB[J]](nnBx[J],nnBy[J])
        if !(CT == Nothing)
            tmpf = intercept[sgiF[J]](tmpf,nnFx[J],nnFy[J],t)
            tmpb = intercept[sgiB[J]](tmpb,nnBx[J],nnBy[J],t)
        end
        
        w[J] = (tmpf + tmpb)/2
    end
    τglobal[I] = 0.1 * (maximum(abs.(u - w))/ maximum(abs.(w)))^2.0
    if isnan(τglobal[I])
        τglobal[I] = zero(TT)
    end
    # isnan(τglobal[I]) ? τglobal[I] = 0.0 : nothing

end


"""
Compute the globalw for parallel blocks where a package is used for bicubic interpolation
"""
function computeglobalw!(u::AbstractArray{TT},uglobal::Vector{Matrix{TT}},τglobal::Vector{TT},Δt::TT,P::Vector{ParallelData{TT,DIM,ParallelGrid{TT,DIM,PMT,AT},GT,BT,IT}},grid::Grid2D{TT,MET},I) where {TT,MET,DIM,AT,GT,BT,IT, PMT<:ParGrid}
    
    Ipt = [BicubicInterpolator(P[I].gridx,P[I].gridy,uglobal[I]) for I in eachindex(uglobal)]

    w = P[I].w

    sgiF = P[I].PGrid.Fplane.subgrid
    sgiB = P[I].PGrid.Bplane.subgrid

    nnFx = P[I].PGrid.Fplane.x
    nnFy = P[I].PGrid.Fplane.y
    nnBx = P[I].PGrid.Bplane.x
    nnBy = P[I].PGrid.Bplane.y

    for J in eachindex(w)
        w[J] = Ipt[sgiF[J]](nnFx[J],nnFy[J])
        w[J] += Ipt[sgiB[J]](nnBx[J],nnBy[J])
        w[J] = w[J]/2
    end

    τglobal[I] = P[I].τ * 0.1 * (maximum(abs.(u - w))/ maximum(abs.(w)))^2.0

end
"""
Compute the globalw for parallel blocks where a custom interpolation scheme is used
"""
function computeglobalw!(u::AbstractArray{TT},uglobal::Vector{Matrix{TT}},τglobal::Vector{TT},Δt::TT,P::Vector{ParallelData{TT,DIM,ParallelGrid{TT,DIM,PMT,AT},GT,BT,IT}},grid::Grid2D{TT,MET},I) where {TT,MET,DIM,AT,GT,BT,IT, PMT<:ParGridLinear{TT,AT,METHOD}} where METHOD
    
    w = P[I].w

    sgiF = P[I].PGrid.Fplane.subgrid
    sgiB = P[I].PGrid.Bplane.subgrid


    nnFx = P[I].PGrid.Fplane.i
    nnFy = P[I].PGrid.Fplane.j
    nnBx = P[I].PGrid.Bplane.i
    nnBy = P[I].PGrid.Bplane.j

    w11F = P[I].PGrid.Fplane.weight11
    w12F = P[I].PGrid.Fplane.weight12
    w21F = P[I].PGrid.Fplane.weight21
    w22F = P[I].PGrid.Fplane.weight22

    w11B = P[I].PGrid.Bplane.weight11
    w12B = P[I].PGrid.Bplane.weight12
    w21B = P[I].PGrid.Bplane.weight21
    w22B = P[I].PGrid.Bplane.weight22

    if METHOD == :NearestNeighbour
        for J in eachindex(u)
            w[J] = uglobal[sgiB[J]][nnBx[J]]
            w[J] += uglobal[sgiF[J]][nnFx[J]]
            w[J] = w[J]/2
        end
    elseif METHOD == :CHS
        error("This should have entered a different `computeglobalw!`. Please check the `typeof(ParallelGrid)`.")
    else # Linear Interpolation
        for J in eachindex(w)
            w[J] = _interpolation(uglobal[sgiB[J]], w11B, w12B, w21B, w22B, nnBx[J], nnBy[J], J)
            w[J] += _interpolation(uglobal[sgiF[J]],w11F, w12F, w21F, w22F, nnFx[J], nnFy[J], J)
            w[J] = w[J]/2
        end
    end
    
    τglobal[I] = P[I].τ * 0.1 * (maximum(abs.(u - w))/ maximum(abs.(w)))^2.0
    if isnan(τglobal[I])
        τglobal[I] = zero(TT)
    end

end





"""
    _interpolation
"""
function _interpolation(u::AT,w11::AT,w12::AT,w21::AT,w22::AT,i::Int,j::Int,I::Int) where {AT}
    w = (u[i,j] * w11[I])
    w += (u[i+1,j] * w12[I])
    w += (u[i,j+1] * w21[I])
    w += (u[i+1,j+1] * w22[I])
    return w
end