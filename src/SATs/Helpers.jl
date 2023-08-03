


function applySAT!(Boundary::Nothing,Block,Prob,mode) end
function applySAT!(SAT::SimultanousApproximationTerm,cache::AT,u::AT,k::AT,mode::SATMode) where AT
    SAT(cache,u,k,mode)
end

function CoupleBlock(BlockA,BlockB,Joint)
end

# Base.eachindex(SB::SATBoundaries{L}) where L = L