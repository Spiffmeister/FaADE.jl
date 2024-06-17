

"""
    Computes the x-y coordinate a point in a domain bounded by `cbottom` `ctop` `cleft` and `cright`.
"""
function coordinate(cbottom::Function,cleft::Function,cright::Function,ctop::Function,u::TT,v::TT) where TT
    S = (one(TT)-v)*cbottom(u) + v*ctop(u) + (one(TT)-u)*cleft(v) + u*cright(v) - 
        (u*v*ctop(one(TT)) + u*(one(TT)-v)*cbottom(one(TT)) + v*(one(TT)-u)*ctop(zero(TT)) + (one(TT)-u)*(one(TT)-v)*cbottom(zero(TT)))
    # if u == 1
    #     @show S
    # end
    return S
end






"""
    meshgrid(S,TT,nx,ny)
Generates the 2D grid of `nx` and `ny` points in a domain given a set of functions which bound a domain.
"""
function meshgrid(TT,cbottom::Function,cleft::Function,cright::Function,ctop::Function,nx::Int,ny::Int)
    # TT = Float64

    u = LinRange(TT(0),TT(1),nx)
    v = LinRange(TT(0),TT(1),ny)

    S(u,v) = coordinate(cbottom,cleft,cright,ctop,u,v)

    
    
    X = zeros(nx,ny)
    Y = zeros(nx,ny)

    for j = 1:ny
        for i = 1:nx
            # println(S(u[i],v[j])," ",u[i],",",v[j])
            coord = S(u[i],v[j])
            X[i,j] = coord[1]
            Y[i,j] = coord[2]
        end
    end
    return X,Y
end
"""
    meshgrid(cbottom::Function,cleft::Function,cright::Function,ctop::Function,nx::Int,ny::Int)
Construct a grid of `nx` by `ny` points in a domain bounded by the functions `cbottom`, `cleft`, `cright`, and `ctop`.
"""
meshgrid(cbottom::Function,cleft::Function,cright::Function,ctop::Function,nx::Int,ny::Int) = meshgrid(Float64,cbottom,cleft,cright,ctop,nx,ny)
"""
    meshgrid(𝒟x::Vector{TT},𝒟y::Vector{TT}) where TT
Generate matrix of coordinates from vectors of coordinates this is useful for packed grids.
"""
function meshgrid(𝒟x::Vector{TT},𝒟y::Vector{TT}) where TT
    nx = length(𝒟x)
    ny = length(𝒟y)
    X = zeros(nx,ny)
    Y = zeros(nx,ny)
    for j = 1:ny
        for i = 1:nx
            X[i,j] = 𝒟x[i]
            Y[i,j] = 𝒟y[j]
        end
    end
    return X,Y
end
"""
    meshgrid(cinner::Function,couter::Function,nx,ny)
Meshgrid for annular domains where the inner and outer boundaries are parameterised boundaries
"""
function meshgrid(cbottom::Function,cleft::Function,nx::Int,ny::Int)
    S(u,v) = coordinate(cbottom,cleft,cright,ctop,u,v)
    
    X = zeros(nx,ny)
    Y = zeros(nx,ny)

    for j = 1:ny
        for i = 1:nx
            X[i,j] = S(u[i],v[j])[1]
            Y[i,j] = S(u[i],v[j])[2]
        end
    end
    return X,Y
end
"""
    meshgrid(inner::Torus,outer::Torus,ζ,nr,nθ)
Take two tori and generate a meshgrid between them at a given toroidal angle `ζ` with `nr` radial points and `nθ` poloidal points.

See also: [`Torus`](@ref)
"""
function meshgrid(inner::Torus,outer::Torus,ζ,nr,nθ)

    AL(u) = inner(0.0,0.0) + u*(outer(0.0,0.0) - inner(0.0,0.0))
    AR(u) = inner(2π,0.0) + u*(outer(2π,0.0) - inner(2π,0.0))

    # cbottom(v) = [0.0,2π*v]
    # ctop(v) = [0.0,2π*v]

    X,Y = meshgrid(AL, u->inner(2π*u,0.0), u->outer(2π*u,0.0), AL, nr, nθ)
    # X,Y = meshgrid(cbottom, u->inner(2π*u,0.0), u->outer(2π*u,0.0), ctop, nr, nθ)
    
    return X,Y
end

