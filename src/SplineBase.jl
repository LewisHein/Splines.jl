import Base: extrema

min_knot_dist_factor = 10000 #knots must be separated by at leas 10000*eps(T) to be considered unique.

function forwardDifference{T<:Number}(x::AbstractArray{T, 1})
    length = size(x, 1)
    return [x[i]-x[i-1] for i in 2:length]
end

type Spline{T<:Number}
    knots::Array{T, 1}
    values::Array{T, 1}
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

	todelete = Array{Int, 1}(0)
	for knot in 2:nKnots
		if knots[knot] == knots[knot-1]
			throw(ArgumentError("Knots must all be unique"))
		end
		if knots[knot]-knots[knot-1] < min_knot_dist_factor*eps(T) #then replace the two knots & values by their averages
			knots[knot] = (knots[knot]+knots[knot-1])/2
			values[knot] = (values[knot]+values[knot-1])/2
			push!(todelete, knot-1)
		end
	end
	
	reverse!(todelete) #reverse the list of knots to delete so that removing one does not affect subsequent removals
	for i in todelete
	    deleteat!(knots, i)
	    deleteat!(values, i)
	end
	nKnots = length(knots)

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
	new(knots, values, a[1:nKnots-1], b, c[1:nKnots-1], d)
    end
end

"""Recalculate a spline with new values, assuming that the corresponding knots are already in s.knots"""
function recalculate!{T}(s::Spline{T}, values::Array{T, 1})
  @assert issorted(s.knots) "Knots are not sorted. This is probably a bug in the Splines library"
	nKnots = length(s.knots)
	nvalues = length(values)

	if nvalues != nKnots
		throw(ArgumentError("Values is of wrong size"))
	end

	if any(isnan(values))
		throw(ArgumentError("values cannot contain NaN"))
	end

	h = forwardDifference(s.knots)
	todelete = Array{Int, 1}(0)
	if minimum(h) < min_knot_dist_factor*eps(T)
		for (knot, knotdiff) in enumerate(h)
			if knotdiff == zero(T)
				throw(ArgumentError("s.knots must all be unique"))
			end
			if knotdiff < min_knot_dist_factor*eps(T) #then replace the two s.knots & values by their averages
				s.knots[knot+1] = (s.knots[knot+1]+s.knots[knot])/2
				values[knot+1] = (values[knot+1]+values[knot])/2
				push!(todelete, knot)
			end
		end
	end
	reverse!(todelete) #reverse the list of knots to delete so that removing one does not affect subsequent removals
	
	for i in todelete
		deleteat!(s.knots, i)
		deleteat!(values, i)
	end
	nKnots = length(s.knots)

	a = values
	if length(todelete) > 0
		h = forwardDifference(s.knots)
	end
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
	s.values = values
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


	bracket = searchsortedlast(s.knots, x)
	if bracket == length(s.knots)
		bracket -= 1
	end

	h = x-s.knots[bracket]

	return @evalpoly(h, s.a[bracket], s.b[bracket], s.c[bracket], s.d[bracket])::T
end

"""calling a spline _f_ with another spline _g_ as argument produces a spline approximation to f(g)"""
function call{T}(f::Spline{T}, g::Spline{T})
  gmax, gmin = extrema(g)
  if !indomain(f, gmax) || !indomain(f, gmin)
    throw(DomainError())
  end

  # approixmating f(g) is best done by first approximating g with first-order polynomials, and then strechting or shrinking intervals of f appropriately
  g_discrete = adaptive_discretize(g)

  #However, it is very necessary that the knots/values in f be preserved in the transformation;
  #otherwise, the transformed spline might bear little resemblance to what it should.
  #The simplest way to do this is, for each knot $f_k$ in f to put a knot $g_x$ into g's list of knots such that g($g_k$) = $f_k$

  newknots = g_discrete[:x]
  newvalues = [f(x) for x in g_discrete[:y]]
  f_warped = Spline(newknots, newvalues)

  #We now have some approximation to f(g(t)), but it probably isn't very good.
  # So we begin going through and trying values halfway between the knots (from g_discrete[:x])
  # If these values differ too greatly from the true f(g(t)), then we add a knot there.
  tol = .001 #FIXME: This is purely arbitrary and probably a bad value
  #FIXME: This loop implicitly assumes that the maximum error occurs halfway between knots.
  #This assumption is probably not too terrible, but it would be better to find the true location
  #of the max error and put the knot there
  maxϵ = typemax(T)
  ϵ = zero(T)
  while maxϵ > tol
    knots_toInsert = Array{Float64, 1}(0)
    values_toInsert = Array{Float64, 1}(0)
    maxϵ = typemin(T)
    for i in 2:length(f_warped.knots)
      x = (f_warped.knots[i]+f_warped.knots[i-1])/2
      ϵ = abs(f_warped(x)-f(g(x)))
      if ϵ > tol
	push!(knots_toInsert, x)
	push!(values_toInsert, f(g(x)))
      end
    end
    insert!(f_warped, knots_toInsert, values_toInsert)
    
    for i in 2:length(f_warped.knots)
      x = (f_warped.knots[i]+f_warped.knots[i-1])/2
      ϵ = abs(f_warped(x)-f(g(x)))
      if ϵ > maxϵ
      	maxϵ = ϵ
      end
    end
  end


  return f_warped
end

"""Add a knot to a spline"""
function insert!{T}(s::Spline{T}, knot::T, value::T)
	knots = s.knots
	values = [s(i)::T for i in knots] #FIXME: why allocate a whole new array? all but one of the values are already stored in s.a
	if knot in knots
		values[findin(knots, knot)[1]] == value #FIXME: this should be assignment
		recalculate!(s, values)
		return nothing
	end

	push!(knots, knot)
	push!(values, value)

	perm = sortperm(knots)
	s.knots = knots[perm]
	values = values[perm]

	recalculate!(s, values)
	return nothing
end

#"""Add a bunch of knots to a spline all at once"""
function insert!{T}(s::Spline{T}, knots::AbstractArray{T, 1}, values::AbstractArray{T, 1})
  if length(knots) != length(values)
    throw(ArgumentError("knots and values must be the same size"))
  end
  if unique(knots) != knots
    throw(ArgumentError("knots must all be unique"))
  end

  newvalues = copy(s.values)

  #If there are knots in _knots_ that are already in the spline ,then just change the value and add
  #that index to a list of knots and values to be deleted
  todelete = Array{Int, 1}(0)
  for (i, knot) in enumerate(knots)
    if knot in s.knots
      duplicate_pos = findin(s.knots, knot)[1]
      newvalues[duplicate_pos] = values[i] #change the value
      push!(todelete, i) #and add the index to a list of indecies to delete
    end
  end
  #delete whatever was in the list
  reverse!(todelete) #we have to delete starting at the end so that a deletion at one index does not throw off future deletions
  for i in todelete
    deleteat!(knots, i)
    deleteat!(values, i)
  end

  append!(s.knots, knots)
  append!(newvalues, values)

  #Now the knots must be sorted to be in increasing order and the values must re-arranged in the same manner
  perm = sortperm(s.knots)
  s.knots = s.knots[perm]
  newvalues = newvalues[perm]

  recalculate!(s, newvalues)
  return nothing
end

"""Move the knot at knotIndex to the location newlocation"""
function moveknot!{T}(s::Spline{T}, knotIndex::Int, newlocation::T)
  if knotIndex > length(s.knots) || knotIndex <1
    throw(BoundsError(s, knotIndex))
  end

  if !indomain(s, newlocation)
    throw(DomainError())
  end

  if newlocation in s.knots
    throw(ArgumentError("newlocation cannot be a pre-existing knot"))
  end

  s.knots[knotIndex] = newlocation

  recalculate!(s, s.values)
end

"""Move the knots at _knotIndecies_ to the new locations in _newlocations_"""
function moveknot!{T}(s::Spline{T}, knotIndecies::AbstractArray{Int, 1}, newlocations::AbstractArray{T, 1})
  if length(knotIndecies) != length(newlocations)
    throw(ArgumentError("knotIndecies and newlocations must have the same length, not $(length(knotIndecies)) and $(length(newlocations))"))
  end

  for index in knotIndecies
    if index > length(s.knots) || index < 1
      throw(BoundsError(s, index))
    end
  end

  for newloc in newlocations
    if !indomain(s, newloc)
      throw(DomainError)
    end
  end

  for (i, index) in enumerate(knotIndecies)
    s.knots[index] = newlocations[i]
  end

  recalculate!(s, s.values)
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

include("cubic_helpers.jl")

function maximum{T}(s::Spline{T})
	maxS = -Inf

	for i in 1:length(s.knots)-1
		maxS = max(maxS, maxcubic(s.knots[i+1]-s.knots[i], s.a[i], s.b[i], s.c[i], s.d[i])[1])
	end

	return maxS
end

function minimum{T}(s::Spline{T})
	minS = typemax(T)

	for i in 1:length(s.knots)-1
		minS = min(minS, mincubic(s.knots[i+1]-s.knots[i], s.a[i], s.b[i], s.c[i], s.d[i])[1])
	end

	return minS
end

function extrema{T}(s::Spline{T})
	minS = typemax(T)
	maxS = typemin(T)

	for i in 1:length(s.knots)-1
		minS = min(minS, mincubic(s.knots[i+1]-s.knots[i], s.a[i], s.b[i], s.c[i], s.d[i])[1])
		maxS = max(maxS, maxcubic(s.knots[i+1]-s.knots[i], s.a[i], s.b[i], s.c[i], s.d[i])[1])
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


"""Find the appropriate step size for discretizing a spline with a uniform mesh maximum absolute error _tol_ due to linear interpolation. If negative, this will be interpreted as a relative error from the maximum value"""
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
function uniform_discretize{T}(s::Spline{T}, tol::Number=-1) #default .001 absolute error tolerated
	mesh = uniform_discretize_mesh(s, tol)

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

"""A helper function to find the maximum absolute error due to linear interpolation of
a cubic spline over a mesh _x_

Parameters:
 - *s* is the spline that is being discretized
 - *x* holds x-values that need to be used between the current knot (which is not easily visible to this function) and the one after it to get a reasonable approximation.
 - *a*, *b*, *c*, and *d* are the polynomial coefficients of the spline at the knot where this is being done.
"""
function find_maxϵ{T}(s::Spline{T}, x::Array{T, 1}, a::T, b::T, c::T, d::T)
  if any(isnan(x))
    throw(ArgumentError("NAN encountered in find_maxϵ"))
  end

  maxϵ = maxϵ_pos = maxϵ_index = typemin(T) #maxcubic(x[1], x[2], 0, s.b[i]-((s(s.knots[i]+x[2])-s(s.knots[i]+x[1]))/(x[2]-x[1])), s.c[i], s.d[i])
	for (j, xj) in take(enumerate(x), length(x)-1)
    lin_approx_slope = (s(x[j+1])-s(xj))/(x[j+1]-xj)
	  maxϵ_new, maxϵ_pos_new = maxabscubic(xj-x[1], x[j+1]-x[1], a-(s(xj)-lin_approx_slope*(xj-x[1])), b-lin_approx_slope, c, d)
    @assert maxϵ_pos != xj-x[1] "$xj $(a-(s(xj)-lin_approx_slope*(xj-x[1])))"

	  if maxϵ_new > maxϵ
	    maxϵ = maxϵ_new
		  maxϵ_pos = maxϵ_pos_new
		  maxϵ_index = j
    end
  end

  return maxϵ, maxϵ_pos, maxϵ_index
end



"""
A helper function for adaptive_discretize_mesh

Parameters:
 - *s* is the spline that is being discretized.
 - *x* holds x-values that need to be used between the current knot (which is not easily visible to this function) and the one after it to get a reasonable approximation.
 - *absTol* is the maximum absolute error that is tolerated
 - *a*, *b*, *c*, and *d* are the polynomial coefficients of the spline at the knot where this is being done.

Return:
 - The new error after adding zero or more values to the mesh *x*

Given _s_ and _x_, this function finds the maximum error associated with using linear
interpolation over the mesh _x_. If this error is too large, then it adds a new value
to _x_ where the max error was and returns the new error.
"""
function densify_mesh!{T}(s::Spline{T}, x::Array{T, 1}, absTol::T, a::T, b::T, c::T, d::T)
  maxϵ, maxϵ_pos, maxϵ_index = find_maxϵ(s, x, a, b, c, d)
  @assert !isnan(maxϵ_pos)
  if maxϵ > absTol
    if maxϵ_pos == x[maxϵ_index+1]
      insert!(x, maxϵ_index+1, (x[maxϵ_index+1]+x[maxϵ_index])/2)
    elseif maxϵ_pos in x
      error("duplicate value $maxϵ_pos x=$x")
    else
      insert!(x, maxϵ_index+1, x[1]+maxϵ_pos)
    end
    #but maxϵ has now changed, so we re-compute it
    maxϵ = find_maxϵ(s, x, a, b, c, d)[1]
  end

  return maxϵ
end

"""
Create an adaptive mesh of points for discretizing a spline _s_.

Parameters:
 - *s* A spline to discretize
 - *tol* The tolerance. If this is negative, then it is assumed to be a relative tolerance (with a few modifications to get around the problem of the spline being 0 somewhere)
 otherwise, it is assumed to be an absolute tolerance.
"""
function adaptive_discretize_mesh{T}(s::Spline{T}, tol::Number=-.01)

	if tol > 0
		absTol = tol
	end
	nKnots = length(s.knots)

	x = Array{Array{T, 1}, 1}(nKnots)

	for i in 1:nKnots-1
		x[i] = [s.knots[i], s.knots[i+1]]
		maxϵ = maxabscubic(0.0, x[i][2]-x[i][1], zero(T), s.b[i]-((s.values[i+1]-s.values[i])/(s.knots[i+1]-s.knots[i])), s.c[i], s.d[i])[1]

		if tol < 0
			absTol = abs(tol)*maxabscubic(zero(T), s.knots[i+1]-s.knots[i], s.a[i], s.b[i], s.c[i], s.d[i])[1]
		end

		while maxϵ > absTol
      maxϵ = densify_mesh!(s, x[i], absTol, s.a[i], s.b[i], s.c[i], s.d[i])
		end
		deleteat!(x[i], length(x[i]))
	end
	x[nKnots] = [s.knots[nKnots]]

	return vcat(x...)::Array{T, 1}
end

function adaptive_discretize{T}(s::Spline{T}, tol::Number=-.01)
	mesh = adaptive_discretize_mesh(s, tol)
	y = [s(x) for x in mesh]

	return Dict(:x => mesh, :y => y)
end


export call
export insert!
export moveknot!
export deleteknotat!
export trimstart!
export trimend!
export trim!
export maximum
export minimum
export extrema
export xmax
export maxabs
export Spline
export discretize
export uniform_discretize
export adaptive_discretize
export indomain
export domain
