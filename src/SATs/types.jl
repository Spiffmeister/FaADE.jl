

"""
    SATMode

Used when the conjugate gradient solver is being used so the solver knows which part of the SAT to call, since ``F`` in ``Ab=F`` contains all data (the right hand side of ``u(0)=f(t)``, but not the solution ``u``). 
"""
struct SATMode{T} end
const DataMode = SATMode{:DataMode}()
const SolutionMode = SATMode{:SolutionMode}()

