__precompile__()
module Splines
import Base: LinAlg.LAPACK.gtsv!, maximum, minimum, maxabs, insert!
include("SplineBase.jl")
include("SplineMath.jl")
include("SplineCalculus.jl")
end # module
