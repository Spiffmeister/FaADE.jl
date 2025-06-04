

"""
    BDF2(DBlock::MultiDataBlock, Δt, t)
Computes the solution to the BDF2 integrator for the problem stored in `DBlock` with timestep `Δt` and time `t`.

Assumes `u_0` and `u_1` have already been computed, assume `b` has the RHS stored
"""
function BDF2(DBlock::MultiDataBlock, Δt, t)

    # Compute the RHS
    # RHS!(:uₙ₊₁.DBlock)

    @. DBlock[1].cache = DBlock[1].uₙ + Δt * DBlock[1].b

    conj_grad!(DBlock,BDF2_Residual())

    if DBlock.SC.converged
        # if converged shuffle the solutions along
        setValue(:u,:uₙ₊₁,DBlock)
        setValue(:uₙ₊₁,:cache,DBlock)
    end
end


"""
    BDF2_Residual(DB::LocalDataBlock{TT},Δt::TT,t::TT)

    Compute the residual for the BDF2 integrator storing the result in rₖ and where `cache = uₙ₊₁`
"""
function BDF2_Residual(DB::LocalDataBlock{TT},Δt::TT,t::TT)

    @. D.rₖ = TT(3) * D.cache - TT(4)*D.u + TT(2)*Δt*D.cache

end