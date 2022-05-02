using MathOptInterface
const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIU = MOI.Utilities
const MOIB = MOI.Bridges

using Test

# First pure backend test
include("MOI/backend-tests.jl")

# Load configs only once
include("MOI/configs.jl")

solverstacks = get_solverstacks()
@testset "Tests with solver backend $stack" for stack in solverstacks
  set_solverstack(stack)
  include("MOI/solvers-$stack.jl")
  include("MOI/wrapper.jl")

  if stack == "GAMS"
    include("export_gms.jl")
    export_gms(mktempdir())

    include("MOI/minlptests.jl")

    # GAMS is really fragile on macos and windows
    if Sys.islinux()
      GC.gc()

      include("MOI/solvers_fragile.jl")
      include("MOI/wrapper.jl")
    end
  end
end
