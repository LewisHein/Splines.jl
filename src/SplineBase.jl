function forwardDifference{T<:Number}(x::AbstractArray{T, 1})
    length = size(x, 1)
    return [x[i]-x[i-1] for i in 2:length]
end

type Spline{T<:Number}
    knots::Array{T, 1}
    a::Array{T, 1}
    b::Array{T, 1}
    c::Array{T, 1}
    d::Array{T, 1}


    function Spline(knots::AbstractArray{T, 1}, values::AbstractArray{T, 1})
	nKnots = length(knots)
    	if nKnots <= 1
		throw(ArgumentError("There must be at least two knots"))
	end

	if length(values) != nKnots
		throw(ArgumentError("Values is of wrong size"))
	end

	if any(isnan(knots))
		throw(ArgumentError("Knots cannot be NaN"))
	end

	if any(isnan(values))
		throw(ArgumentError("Values cannot be NaN"))
	end

        if !issorted(knots)
		throw(ArgumentError("Knots must be sorted"))
	end

	for knot in 2:nKnots
		if knots[knot] == knots[knot-1]
			throw(ArgumentError("Knots must all be unique"))
		end
	end

	a = values
	h = forwardDifference(knots)
	df = forwardDifference(values) ./ h
# Note: due the way that the gtsv! solver works, when solving A*X = B, this matrix can be X and then be overwritten By B. Therefore, the name of c is logical
	c = Array{Float64, 1}([0])
	append!(c, 3*forwardDifference(df))
	append!(c, [0])

	lowerDiag::Array{Float64, 1} = h[1:nKnots-2]
	append!(lowerDiag, [0.0]);

	diag = Array{Float64}([1])
	append!(diag, [2*(h[i]+h[i-1]) for i in 2:nKnots-1])
	append!(diag, [1.0])

	upperDiag = Array{Float64, 1}([0.0])
	append!(upperDiag, h[2:nKnots-1])

	#LAPACK.gtsv!([h[1:nKnots-2]; 0], [1; [2*(h[i]+h[i-1]) for i in 2:nKnots-1]; 1], [0; h[2:nKnots-1]], c;)
	LAPACK.gtsv!(lowerDiag, diag, upperDiag, c)
	d = [(c[i+1]-c[i])/(3*h[i]) for i in 1:nKnots-1]
	b = [df[i]-((h[i]/3)*((2*c[i])+c[i+1])) for i in 1:nKnots-1]
	new(knots, a[1:nKnots-1], b, c[1:nKnots-1], d)
    end
end

"""Recalculate a spline with new values"""
function recalculate!{T}(s::Spline{T}, values::Array{T, 1})
	nKnots = length(s.knots)
	nvalues = length(values)

	if nvalues != nKnots
		throw(ArgumentError("Values is of wrong size"))
	end

	if any(isnan(values))
		throw(ArgumentError("values cannot contain NaN"))
	end

	a = values
	h = forwardDifference(s.knots)
	df = forwardDifference(values) ./ h
# Note: due the way that the gtsv! solver works, when solving A*X = B, this matrix can be X and then be overwritten By B. Therefore, the name of c is logical
	c = Array{Float64, 1}([0])
	append!(c, 3*forwardDifference(df))
	append!(c, [0])

	lowerDiag::Array{Float64, 1} = h[1:nKnots-2]
	append!(lowerDiag, [0.0]);

	diag = Array{Float64}([1])
	append!(diag, [2*(h[i]+h[i-1]) for i in 2:nKnots-1])
	append!(diag, [1.0])

	upperDiag = Array{Float64, 1}([0.0])
	append!(upperDiag, h[2:nKnots-1])

	#LAPACK.gtsv!([h[1:nKnots-2]; 0], [1; [2*(h[i]+h[i-1]) for i in 2:nKnots-1]; 1], [0; h[2:nKnots-1]], c;)
	LAPACK.gtsv!(lowerDiag, diag, upperDiag, c)
	d = [(c[i+1]-c[i])/(3*h[i]) for i in 1:nKnots-1]
	b = [df[i]-((h[i]/3)*((2*c[i])+c[i+1])) for i in 1:nKnots-1]
	s.a = a[1:nKnots-1]
	s.b = b
	s.c = c[1:nKnots-1]
	s.d = d
#	new(knots, a[1:nKnots-1], b, c[1:nKnots-1], d)
	return nothing
end

function indomain{T}(s::Spline{T}, x::T)
	return (x>=s.knots[1] && x<=s.knots[end])
end

function domain{T}(s::Spline{T})
	return (s.knots[1], s.knots[end])
end

Spline{T<:Number}(knots::AbstractArray{T, 1}, values::AbstractArray{T, 1}) = Spline{T}(knots, values)

function call{T}(s::Spline{T}, x::T)
	if !indomain(s, x)
		throw(DomainError())
		#error("x out of range")
	end

	bracket = 0
	for i in size(s.knots, 1)-1:-1:1
		if x >= s.knots[i]
			bracket = i
			break
		end
	end

	h = x-s.knots[bracket]

	return @evalpoly(h, s.a[bracket], s.b[bracket], s.c[bracket], s.d[bracket])
end

"""Add a knot to a spline"""
function insert!{T}(s::Spline{T}, knot::T, value::T)
	knots = s.knots
	values = [s(i)::T for i in knots] #FIXME: why allocate a whole new array? all but one of the values are already stored in s.a
	if knot in knots
		values[findin(knots, knot)[1]] == value
		recalculate!(s, values)
		return nothing
	end

	push!(knots, knot)
	push!(values, value)

	perm = sortperm(knots)
	s.knots = knots[perm]
	values = values[perm]

	recalculate!(s, values)
end

"""Delete a knot or knots from a spline"""
function deleteknotat!{T}(s::Spline{T}, knotIndex::Integer...)
	knots = s.knots
	values = Array{T, 1}(length(knots))
	for (i, x) in enumerate(knots)
		values[i] = s(x)
	end

	deleteat!(knots, knotIndex)
	deleteat!(values, knotIndex)

	recalculate!(s, values)
end

"""Trim the beginning (everything less than cutoff) off of a spline"""
function trimstart!{T}(s::Spline{T}, cutoff::T)
	if !indomain(s, cutoff)
		throw(DomainError())
	end

	newknot = cutoff
	newvalue = s(newknot)
	insert!(s, newknot, newvalue)

	indeciesToTrim = find(x->x < newknot, s.knots)
	deleteknotat!(s, indeciesToTrim...)
end

"""Trim the end (everything greater than cutoff) off of a spline"""
function trimend!{T}(s::Spline{T}, cutoff::T)
	if !indomain(s, cutoff)
		throw(DomainError())
	end

	newknot = cutoff
	newvalue = s(newknot)
	insert!(s, newknot, newvalue)

	indeciesToTrim = find(x->x > newknot, s.knots)
	deleteknotat!(s, indeciesToTrim...)
end

"""Trim the start and the end (everything less than startx or greater than endx) off of a spline"""
function trim!{T}(s::Spline{T}, startx::T, endx::T)
	insert!(s, startx, s(startx))
	insert!(s, endx, s(endx))
	indeciesToTrim = find(x->x<startx || x>endx, s.knots)
	deleteknotat!(s, indeciesToTrim...)
end


#FIXME: NaN if d = 0
function maxcubic(h, a, b, c, d) #a+bx+cx^2+dx^3
	#We use the quadratic formula on the derivative
	discriminant = (2c)^2-4*(3d*b)
	if discriminant < 0
		return -Inf
	end

	root1 = (-2*c-sqrt(discriminant))/(6*d)
	root2 = (-2*c+sqrt(discriminant))/(6*d)

	val1 = @evalpoly(0, a, b, c, d)
	val2 = @evalpoly(root1, a, b, c, d)
	val3 = @evalpoly(root2, a, b, c, d)
	val4 = @evalpoly(h, a, b, c, d)

	if root1 < 0 || root1 >h
		val2 = -Inf
	end
	if root2 < 0 || root2 > h
		val3 = -Inf
	end

	theMax = max(val1, val2, val3, val4)

	maxPos = -Inf
	if val1 == theMax
		maxPos = 0
	elseif val2 == theMax
		maxPos = root1
	elseif val3 == theMax
		maxPos = root2
	elseif val4 == theMax
		maxPos = h
	else
		error("Failed to find the max in maxcubic. This is probably a bug in Splines.jl")
	end

	return (theMax, maxPos)
end


function mincubic(h, a, b, c, d) #a+bx+cx^2+dx^3
	#We use the quadratic formula on the derivative
	discriminant = (2*c)^2-4*(3*d*b)
	if discriminant < 0
		return Inf
	end

	root1 = (-2*c-sqrt(discriminant))/(6*d)
	root2 = (-2*c+sqrt(discriminant))/(6*d)

	val1 = @evalpoly(0, a, b, c, d)
	val2 = @evalpoly(root1, a, b, c, d)
	val3 = @evalpoly(root2, a, b, c, d)
	val4 = @evalpoly(h, a, b, c, d)

	if root1 < 0 || root1 > h
		val2 = Inf
	end

	if root2 < 0 || root2 > h
		val3 = Inf
	end

	theMin = min(val1, val2, val3, val4)

	minPos = Inf
	if val1 == theMin
		minPos = 0
	elseif val2 == theMin
		minPos = root1
	elseif val3 == theMin
		minPos = root2
	elseif val4 == theMin
		minPos = h
	else
		error("Failed to find the minium. This is probably a bug in Splines.jl")
	end

	return (theMin, minPos)
end


function maximum{T}(s::Spline{T})
	maxS = -Inf

	for i in 1:length(s.knots)-1
		maxS = max(maxS, maxcubic(s.knots[i+1]-s.knots[i], s.a[i], s.b[i], s.c[i], s.d[i])[1])
	end

	return maxS
end

function minimum{T}(s::Spline{T})
	minS = Inf

	for i in 1:length(s.knots)-1
		minS = min(minS, mincubic(s.knots[i+1]-s.knots[i], s.a[i], s.b[i], s.c[i], s.d[i])[1])
	end

	return minS
end

function extrema{T}(s::Spline{T})
	minS = Inf
	maxS = -Inf

	for i in 1:length(s.knots)-1
		minS = min(minS, mincubic(s.knots[i+1]-s.knots[i], s.a[i], s.b[i], s.c[i], s.d[i]))
		maxS = max(maxS, maxcubic(s.knots[i+1]-s.knots[i], s.a[i], s.b[i], s.c[i], s.d[i]))
	end

	return (minS, maxS)
end


"""Similar to indmax for arrays, return the x-position of the maximum of the given spline"""
function xmax{T}(s::Spline{T})
	theMax = maximum(s)
	maxPos = 0
	curMax = (0, 0)
	for i in 1:length(s.knots)-1
		curMax = maxcubic(s.knots[i+1]-s.knots[i], s.a[i], s.b[i], s.c[i], s.d[i])
		if curMax[1] == theMax
			maxPos = s.knots[i] + curMax[2]
			break
		end
	end

	return maxPos
end


function maxabs(s::Spline)
	return max(abs(maximum(s)), abs(minimum(s)))
end


"""Compute the maximum absolute value of the second derivative of s"""
function maxabs2deriv(s::Spline)
	maxA2D = sum = 0
	for i in 1:length(s.knots)-1
		maxA2D = max(maxA2D, abs(sum += (2*s.c[i]+6*s.d[i]*(s.knots[i+1]-s.knots[i]))))
	end

	return maxA2D
end


function derivative{T}(s::Spline{T})
	nSplines = size(s.a, 1)

	sDeriv = s
	for i in 1:nSplines
		sDeriv.a[i] = sDeriv.b[i]
		sDeriv.b[i] = sDeriv.c[i]*2
		sDeriv.c[i] = sDeriv.d[i]*3
		sDeriv.d[i] = 0
	end

	return sDeriv
end

"""Find the appropriate step size for discretizing a spline with maximum absolute error _tol_ due to linear interpolation. If negative, this will be interpreted as a relative error from the maximum value"""
function find_stepsize{T}(s::Spline{T}, tol::Number)
	if tol < 0
		tol = maximum(s)*abs(tol)
	end
	minValues = 100
	maxa2d = maxabs2deriv(s)
	if maxa2d > 0
		stepSize = .9*min(sqrt(8*tol/maxa2d), (s.knots[end]-s.knots[1])/minValues)
	else
		stepSize = (s.knots[end]-s.knots[1])/minValues
	end

	return stepSize
end

"""
Create an evenly spaced mesh of points to discretize the spline so that the maximum absolute error due to linear interpolation is <= tol.
If this is too small, return a mesh with 100 points. Note, this may use horrendous amounts of memory. Consider using adaptive_discretize_mesh instead
"""
function uniform_discretize_mesh{T}(s::Spline{T}, tol::Number=-1)
	stepSize = find_stepsize(s, tol)
	mesh = s.knots[1]:stepSize:s.knots[end]
	return mesh
end


"""Turn a spline object into a discrete list of values on an automatically created uniform mesh"""
function discretize{T}(s::Spline{T}, tol::Number=-1) #default .001 absolute error tolerated
	mesh = discretize_mesh(s, tol)

	return discretize(s, mesh)
end

"""Turn a spline object into a discrete list of values at the mesh points given by _mesh_"""
function discretize{T}(s::Spline{T}, mesh::AbstractArray{T})
	discretized = Array{T, 1}(length(mesh))
	for (i, x) in enumerate(mesh)
		discretized[i] = s(x)
	end

	return discretized
end


export call
export insert!
export deleteknotat!
export trimstart!
export trimend!
export trim!
export maximum
export minimum
export extrema
export xmax
export maxabs
export maxabs2deriv
export Spline
export derivative
export discretize
export indomain
export domain
