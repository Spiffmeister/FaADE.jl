using SafeTestsets


@safetestset "Curvilinear Derivatives"  begin include("Derivatives/CurvilinearDerivatives.jl")  end
@safetestset "First Derivative"         begin include("Derivatives/FirstDerivative.jl")         end
@safetestset "First Derivative"         begin include("Derivatives/SecondDerivative.jl")        end

@safetestset "Grid Construction"        begin include("Grid/Gridconstruction.jl")               end

@safetestset "Full 1D run"              begin include("runs/1DFullTest.jl")                     end

@safetestset "Full 2D run"              begin include("runs/2DFullTest.jl")                     end


