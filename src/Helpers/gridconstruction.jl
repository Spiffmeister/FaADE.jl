


function coordinate(cbottom::Function,cleft::Function,cright::Function,ctop::Function,u::TT,v::TT) where TT
    S = (one(TT)-v)*cbottom(u) + v*ctop(u) + (one(TT)-u)*cleft(v) + u*cright(v) - 
        (u*v*ctop(one(TT)) + u*(one(TT)-v)*cbottom(one(TT)) + v*(one(TT)-u)*ctop(zero(TT)) + (one(TT)-u)*(one(TT)-v)*cbottom(zero(TT)))
    return S
end




"""
    meshgrid(TT,cbottom::Function,cleft::Function,cright::Function,ctop::Function,nx::Int,ny::Int)
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
            X[i,j] = S(u[i],v[j])[1]
            Y[i,j] = S(u[i],v[j])[2]
        end
    end
    return X,Y
end
"""
    meshgrid(S,nx,ny)
"""
meshgrid(cbottom::Function,cleft::Function,cright::Function,ctop::Function,nx::Int,ny::Int) = meshgrid(Float64,cbottom,cleft,cright,ctop,nx,ny)






"""
    circulargrid(R::Float64,n::Int)
"""
function circulargrid(R::Float64,n::Int)
    
    # each of the square edges
    squpper(u)  = [u/2 - 0.25,  0.25]
    sqleft(v)   = [-0.25,       v/2 - 0.25]
    sqright(v)  = [0.25,        v/2 - 0.25]
    sqbottom(u)    = [u/2 - 0.25,  -0.25]


    # The 45 degree lines
    sq45left(u) = [u/2 - 0.25,  u/2 - 0.25]
    sq45right(u)= [0.25 - u/2,  u/2 - 0.25]
    sq45top(u)  = [u/2 - 0.25,  0.25 - u/2]
    sq45bottom(u) = [u/2 - 0.25,  u/2 - 0.25]



    # Top annulus
    AUbottom(u) = [u/2 - 0.25, 0.25]
    AUleft(v)   = v*[cos(3π/4) + 0.25, sin(3π/4) - 0.25] + [-0.25, 0.25]
    AUright(v)  = v*[cos(π/4) - 0.25, sin(π/4) - 0.25] + [0.25, 0.25]
    AUtop(u)    = [cos(u*(π/4 - 3π/4) + 3π/4), sin(u*(π/4 - 3π/4) + 3π/4)]

    # Left annulus
    ALbottom(u) = u*[-cos(3π/4) + 0.25, -sin(3π/4) + 0.25] + [-0.25, 0.25]
    ALleft(v)   = [cos(v*(3π/4 - 5π/4) + 5π/4), sin(v*(3π/4 - 5π/4) + 5π/4)]
    ALright(v)  = [-0.25, v/2 - 0.25]
    ALtop(u)    = u*[0.25 + cos(3π/4), 0.25 + sin(π/4)] + [-cos(3π/4), -sin(π/4)]

    # Bottom annulus
    ABbottom(u) = [cos(u*(7π/4 - 5π/4) + 5π/4), sin(u*(7π/4 - 5π/4) + 5π/4)]
    ABleft(v)   = v*[-0.25 - cos(5π/4), -0.25 - sin(5π/4)] + [cos(5π/4), sin(5π/4)]
    ABright(v)  = v*[0.25 - cos(7π/4), -0.25 - sin(7π/4)] + [cos(7π/4), sin(7π/4)]
    ABtop(u)    = [u/2 - 0.25, -0.25]

    # Right annulus
    ARbottom(u) = u*[cos(7π/4) - 0.25, sin(7π/4) + 0.25] + [0.25, -0.25]
    ARleft(v)   = [0.25, v/2 - 0.25]
    ARright(v)  = [cos(v*(9π/4 - 7π/4) + 7π/4), sin(v*(9π/4 - 7π/4) + 7π/4)]
    ARtop(u)    = u*[cos(π/4) - 0.25, sin(π/4) - 0.25] + [0.25, 0.25]



    sq = Grid2D(sqbottom,sqleft,sqright,squpper,nx,ny)
    AU = Grid2D(AUbottom,AUleft,AUright,AUtop,nx,ny)
    AL = Grid2D(ALbottom,ALleft,ALright,ALtop,nx,ny)
    AB = Grid2D(ABbottom,ABleft,ABright,ABtop,nx,ny)
    AR = Grid2D(ARbottom,ARleft,ARright,ARtop,nx,ny)
    
    
    glayout = ([(2,Up),(3,Left),(4,Down),(5,Right)],
                [(1,Down),(3,Left),(5,Right)],
                [(1,Right),(2,Up),(4,Down)],
                [(1,Up),(3,Right),(5,Right)],
                [(1,Left),(2,Left),(4,Down)])
    
    
    G = GridMultiBlock((sq,AU,AL,AB,AR),glayout)
    
    

    return x,y
end



