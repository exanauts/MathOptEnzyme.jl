module MathOptEnzyme

using MathOptInterface
using Enzyme

import Base.Meta: isexpr
import SparseArrays
import Symbolics

const MOI = MathOptInterface

struct EnzymeBackend <:
    MOI.Nonlinear.AbstractAutomaticDifferentiation end

include("nonlinear_oracle.jl")

end
