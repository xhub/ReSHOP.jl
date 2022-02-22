MOI.supports(::Optimizer, ::MOI.VariableName, ::Type{MOI.VariableIndex}) = false

MOI.supports_constraint(::Optimizer, ::Type{<:ALLV}, ::Type{<:SS}) = true
MOI.supports_constraint(::Optimizer, ::Type{<:ALLV}, ::Type{<:NCS}) = true

MOI.supports_constraint(::Optimizer, ::Type{<:SF}, ::Type{<:LS}) = true
MOI.supports_constraint(::Optimizer, ::Type{VAF}, ::Type{<:VLS}) = true
MOI.supports_constraint(::Optimizer, ::Type{VOV}, ::Type{<:VLS}) = true
MOI.supports_constraint(::Optimizer, ::Type{VOV}, ::Type{MOI.Complements}) = true
MOI.supports_constraint(::Optimizer, ::Type{VAF}, ::Type{MOI.Complements}) = true

# TODO:
MOI.supports(::Optimizer, ::MOI.ConstraintName, ::Type{MOI.ConstraintIndex}) = false
# TODO we can support those
MOI.supports_constraint(::Optimizer, ::Type{VOV}, ::Type{NPC}) = false

MOI.supports(::Optimizer, ::MOI.VariablePrimalStart, ::Type{MOI.VariableIndex}) = true
MOI.supports(::Optimizer, ::MOI.ConstraintDualStart, ::Type{MOI.ConstraintIndex}) = true

MOI.supports(::Optimizer, ::MOI.NLPBlock) = true

MOI.supports(::Optimizer, ::MOI.ObjectiveFunction{MOI.SingleVariable}) = true
MOI.supports(::Optimizer, ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}) = true
MOI.supports(::Optimizer, ::MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}}) = true
MOI.supports(::Optimizer, ::MOI.ObjectiveSense) = true


