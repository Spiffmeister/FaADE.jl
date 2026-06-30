

function squared_circle(grid_points::Vector{Tuple{Int,Int}}; circle_radius=1.0, square_width=0.5, dilation=0.0, radial_packing=v->v)

    g0b = square_width/2

    D1 = Grid2D(
        u -> [-g0b,-g0b] +  u*([g0b,-g0b] - [-g0b,-g0b]) + u*(1-u)*[0.0,-dilation],
        v -> [-g0b,-g0b] +  v*([-g0b,g0b] - [-g0b,-g0b]) + v*(1-v)*[-dilation,0.0],
        v -> [g0b,-g0b] +   v*([g0b,g0b] - [g0b,-g0b]) + v*(1-v)*[dilation,0.0],
        u -> [-g0b,g0b] +   u*([g0b,g0b] - [-g0b,g0b]) + u*(1-u)*[0.0,dilation],
        grid_points[1]...
    )

    Tor = Torus([circle_radius],[circle_radius],[1],[0])

    # Right domain
    D2 = Grid2D(#u->[0.25, -u*0.5 + 0.25], # Bottom
        u->[g0b,g0b] + u*([g0b,-g0b] - [g0b,g0b]) + u*(1-u)*[dilation,0.0],
        v->v*(Tor(π/4,0.0) - [g0b,g0b]) + [g0b,g0b], # Left
        v->v*(Tor(7π/4,0.0) + [-g0b, g0b]) + [g0b, -g0b], # Right
        u->Tor(u*(7π/4 - 9π/4) + 9π/4,0.0), # Top
        grid_points[2]...,
        stretchv=radial_packing
    )

    # Top domain
    D3 = Grid2D(#u->[u*0.5 - g0b, g0b], # Bottom
        u->[-g0b,g0b] + u*([g0b,g0b] - [-g0b,g0b]) + u*(1-u)*[0.0,dilation],
        v->v*(Tor(3π/4,0.0) + [g0b,-g0b]) + [-g0b,g0b], # Left
        v->v*(Tor(π/4,0.0) - [g0b,g0b]) + [g0b,g0b], # Right
        u->Tor(u*(π/4 - 3π/4) + 3π/4,0.0), # Top
        grid_points[3]...,
        stretchv=radial_packing
    )

    # Left domain
    D4 = Grid2D(#u->[-g0b,u*0.5 - g0b],
        u->[-g0b,-g0b] + u*([-g0b,g0b] - [-g0b,-g0b]) + u*(1-u)*[-dilation,0.0],
        v->v*(Tor(5π/4,0.0) - [-g0b, -g0b]) + [-g0b, -g0b],
        v->v*(Tor(3π/4,0.0) - [-g0b,g0b]) + [-g0b,g0b],
        u->Tor(u*(3π/4 - 5π/4) + 5π/4,0.0),
        grid_points[4]...,
        stretchv=radial_packing
    )

    # Bottom domain
    D5 = Grid2D(#u->[-u*0.5 + g0b, -g0b],
        u->[-g0b,-g0b] + u*([g0b,-g0b] - [-g0b,-g0b]) + u*(1-u)*[0.0,-dilation],
        v->v*(Tor(7π/4,0.0) - [g0b,-g0b]) + [g0b, -g0b],
        v->v*(Tor(5π/4,0.0) - [-g0b,-g0b]) + [-g0b, -g0b],
        u->Tor(u*(5π/4 - 7π/4) + 7π/4, 0.0),
        grid_points[5]...,
        stretchv=radial_packing
    )

    joints = (
        (Joint(2,Right),Joint(3,Up),Joint(4,Left),Joint(5,Down)),
        (Joint(1,Down),Joint(3,Left),Joint(5,Right)),
        (Joint(1,Down),Joint(4,Left),Joint(2,Right)),
        (Joint(1,Down),Joint(5,Left),Joint(3,Right)),
        (Joint(1,Down),Joint(2,Left),Joint(4,Right))
        )

    Dom = GridMultiBlock((D1,D2,D3,D4,D5),joints)

    return Dom
end