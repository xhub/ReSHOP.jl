
# Type of the equation and variable indices
const RHP_IDXT                = Cint

# Model types
const RHP_MDL_GAMS            = Cuint(0)  # GAMS context
const RHP_MDL_RHP             = Cuint(1)  # internal context
const RHP_MDL_JULIA           = Cuint(2)  # Julia context
const RHP_MDL_AMPL            = Cuint(3)  # AMPL context


# Mapping types
const RHP_EQ_UNSET            = Cuint(0)  #  Equation type unset 
const RHP_EQ_MAPPING          = Cuint(1)  #  Mapping (objective fn, functional part of VI, ...)
const RHP_EQ_CONE_INCLUSION   = Cuint(2)  #  Inclusion in a cone (usual constraint) 
const RHP_EQ_BOOLEAN          = Cuint(3)  #  Boolean relation 

# Variable types
const RHP_VARTYPE_X           = Cuint(0)   # continuous variable
const RHP_VARTYPE_B           = Cuint(1)   # binary variable
const RHP_VARTYPE_I           = Cuint(2)   # integer variable
const RHP_VARTYPE_S1          = Cuint(3)   # special order set 1 (SOS1)
const RHP_VARTYPE_S2          = Cuint(4)   # special order set 2 (SOS2)
const RHP_VARTYPE_SC          = Cuint(5)   # semi-continuous variable
const RHP_VARTYPE_SI          = Cuint(6)   # semi-integer variable
const RHP_VARTYPE_IND         = Cuint(7)   # indicator variable
const RHP_VARTYPE_POLYHEDRAL  = Cuint(8)   # variable in a polyhedral cone
const RHP_VARTYPE_SOC         = Cuint(9)   # variable in a SOC
const RHP_VARTYPE_RSOC        = Cuint(10)  # variable in a rotated SOC
const RHP_VARTYPE_EXP         = Cuint(11)  # variable in a EXP cone
const RHP_VARTYPE_DEXP        = Cuint(12)  # variable in a dual EXP
const RHP_VARTYPE_POWER       = Cuint(13)  # variable in a POWER cone
const RHP_VARTYPE_DPOWER      = Cuint(14)  # variable in a dual POWER cone


const RHP_CONE_NONE         = Cuint(0)    # Unset/non-existent */
const RHP_CONE_R_PLUS       = Cuint(1)    # Non-negative real \f$\mathbb{R}_+\f$ */
const RHP_CONE_R_MINUS      = Cuint(2)    # Non-positive real \f$\mathbb{R}_-\f$  */
const RHP_CONE_R            = Cuint(3)    # Real \f$\mathbb{R}\f$ */
const RHP_CONE_0            = Cuint(4)    # Zero singleton */
const RHP_CONE_POLYHEDRAL   = Cuint(5)    # Zero singleton */
const RHP_CONE_SOC          = Cuint(6)    # Second Order cone */
const RHP_CONE_RSOC         = Cuint(7)    # Rotated Second Order cone */
const RHP_CONE_EXP          = Cuint(8)    # Exponential cone */
const RHP_CONE_DEXP         = Cuint(9)    # Dual Exponential cone */
const RHP_CONE_POWER        = Cuint(10)   # Power cone */
const RHP_CONE_DPOWER       = Cuint(11)   # Dual Power cone */

const RHP_Double   =  Cint(0)
const RHP_Integer  =  Cint(1)
const RHP_Boolean  =  Cint(2)
const RHP_String   =  Cint(3)

const RHP_BASIS_STATUS_LOWER      = Cint(0)
const RHP_BASIS_STATUS_UPPER      = Cint(1)
const RHP_BASIS_STATUS_BASIC      = Cint(2)
const RHP_BASIS_STATUS_SUPERBASIC = Cint(3)
const RHP_BASIS_STATUS_INVALID    = Cint(4)

const RHP_INDEX_NA      = -2
const RHP_INDEX_INVALID = -1

const solver_stat = [
   :Optimal,
   :IterationLimit,
   :TimeLimit,
   :SolverInterrupt,
   :EvalError,
   :Capability,
   :LicenseIssue,
   :UserInterrupt,
   :SetupError,
   :SolverError,
   :InternalError,
   :SolveSkipped,
   :SystemError
]

const model_stat = [
   :OptimalGlobal,
   :OptimalLocal,
   :Unbounded,
   :InfeasibleGlobal,
   :InfeasibleLocal,
   :InfeasibleIntermed,
   :Feasible,
   :Integer,
   :NonIntegerIntermed,
   :IntegerInfeasible,
   :LicenseError,
   :ErrorUnknown,
   :ErrorNoSolution,
   :NoSolutionReturned,
   :SolvedUnique,
   :Solved,
   :SolvedSingular,
   :UnboundedNoSolution,
   :InfeasibleNoSolution
]
