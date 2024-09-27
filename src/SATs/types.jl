
"""
    SimultanousApproximationTerm{Type}
"""
abstract type SimultanousApproximationTerm{Type} end



"""
    SATMode{T}

Used with the conjugate gradient solver so the appropriate function call can be used.
    
Values are 
- `DataMode` for applying boundary data,
- `SolutionMode` for applying the term which only includes `u` such as `u-g` in Dirichlet conditions,
- `ExplicitMode` currently no solver implemented uses this.
"""
struct SATMode{T} end
const DataMode = SATMode{:DataMode}()
const SolutionMode = SATMode{:SolutionMode}()
const ExplicitMode = SATMode{:ExplicitMode}()





