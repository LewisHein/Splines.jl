#=
A file to hold the helper functions used for working with cubic polynomials
=#

function maxcubic{T}(x0::T, x1::T, a::T, b::T, c::T, d::T) #a+bx+cx^2+dx^3
  if x0 == x1
    return (x0, @evalpoly(x0, a, b, c, d))
  end
  if d == 0 #then the max is at an endpoint. #FIXME: I don't think so
    endpt1_val = @evalpoly(x0, a, b, c, d)
    endpt2_val = @evalpoly(x1, a, b, c, d)
    return endpt2_val > endpt1_val ? (endpt2_val, x1) : (endpt1_val, x0)
  end

  #Otherwise, we use the quadratic formula on the derivative
	discriminant = (c)^2-(3d*b)
	if discriminant >= 0 #then the max may be  in [x0, x1]
		root1 = (-sqrt(discriminant)-c)/(3*d)
		root2 = (sqrt(discriminant)-c)/(3*d)
	end


	val1 = @evalpoly(x0, a, b, c, d)::T
	if discriminant >= 0
		val2 = @evalpoly(root1, a, b, c, d)::T
		val3 = @evalpoly(root2, a, b, c, d)::T
	else
		val2 = typemin(T)
		val3 = typemin(T)
	end
	val4 = @evalpoly(x1, a, b, c, d)::T

	if discriminant >= 0
		if root1 < x0 || root1 >x1
			val2 = typemin(T)
		end
		if root2 < x0 || root2 > x1
			val3 = typemin(T)
		end
	end

	theMax = max(val1, val2, val3, val4)

	maxPos = typemin(T)
	if val1 == theMax
		maxPos = x0
	elseif val2 == theMax
		maxPos = root1
	elseif val3 == theMax
		maxPos = root2
	elseif val4 == theMax
		maxPos = x1
	else
		println("Failed to find the max in maxcubic. This is probably a bug in Splines.jl")
		return nothing
	end

	return (theMax, maxPos)
end

function maxcubic{T}(x::T, a::T, b::T, c::T, d::T) #search on [0, x]
	return maxcubic(zero(T), x, a, b, c, d)
end


#FIXME: NaN if d = 0
function maxabscubic{T}(x0::T, x1::T, a::T, b::T, c::T, d::T) #a+bx+cx^2+dx^3
  @assert !isnan(x0)
  @assert !isnan(x1)
  @assert !isnan(a)
  @assert !isnan(b)
  @assert !isnan(c)
  @assert !isnan(d)
  #print("x0, x1, a, b, c, d: ", x0, a, b, c, d)
  #print("\n")
  if x0 == x1 #then there is no real max, but we do the sanest possible thing
    return (x0, @evalpoly(x0, a, b, c, d))
  end
  if d == 0 #then the max is at an endpoint. #FIXME: I don't think so
    endpt1_val = @evalpoly(x0, a, b, c, d)
    endpt2_val = @evalpoly(x1, a, b, c, d)
    return endpt2_val > endpt1_val ? (endpt2_val, x1) : (endpt1_val, x0)
  end

	#otherwise, we use the quadratic formula on the derivative
	discriminant = (c)^2-(3d*b)
	if discriminant >= 0 #then the max may be  in [x0, x1]
		root1 = (-sqrt(discriminant)-c)/(3*d)
		root2 = (sqrt(discriminant)-c)/(3*d)
	end


	val1 = abs(@evalpoly(x0, a, b, c, d)::T)
	if discriminant >= 0
		val2 = abs(@evalpoly(root1, a, b, c, d)::T)
		val3 = abs(@evalpoly(root2, a, b, c, d)::T)
	else
		val2 = typemin(T)
		val3 = typemin(T)
	end
	val4 = abs(@evalpoly(x1, a, b, c, d)::T)

	if discriminant >= 0
		if root1 < x0 || root1 >x1
			val2 = typemin(T)
		end
		if root2 < x0 || root2 > x1
			val3 = typemin(T)
		end
	end

	theMax = max(val1, val2, val3, val4)

  #	println(val1, val2, val3, val4, theMax)
	maxPos = typemin(T)
	if val1 == theMax
		maxPos = x0
	elseif val2 == theMax
		maxPos = root1
	elseif val3 == theMax
		maxPos = root2
	elseif val4 == theMax
		maxPos = x1
  #=	else
    println("Failed to find the max in maxabscubic. This is probably a bug in Splines.jl")
		return nothing=#
	end

	return (theMax, convert(T, maxPos))
end

function maxabscubic{T}(x::T, a::T, b::T, c::T, d::T) #search on [0, x]
	return maxcubic(zero(T), x, a, b, c, d)
end

function mincubic{T}(x0::T, x1::T, a::T, b::T, c::T, d::T) #a+bx+cx^2+dx^3
  if x0 == x1
    return (x0, @evalpoly(x0, a, b, c, d))
  end
  if d == 0 #then the min is at an endpoint. #FIXME: I don't think so
    endpt1_val = @evalpoly(x0, a, b, c, d)
    endpt2_val = @evalpoly(x1, a, b, c, d)
    return endpt2_val < endpt1_val ? (endpt2_val, x1) : (endpt1_val, x0)
  end

	#otherwise, we use the quadratic formula on the derivative
	discriminant = (c)^2-(3*d*b)
	if discriminant >= 0
		root1 = (-sqrt(discriminant)-c)/(3*d)
		root2 = (sqrt(discriminant)-c)/(3*d)
	end


	val1 = @evalpoly(x0, a, b, c, d)
	if discriminant >= 0
		val2 = @evalpoly(root1, a, b, c, d)
		val3 = @evalpoly(root2, a, b, c, d)
	else
		val2 = typemax(T)
		val3 = typemax(T)
	end
	val4 = @evalpoly(x1, a, b, c, d)

	if discriminant >= 0
		if root1 < x0 || root1 > x1
			val2 = typemax(T)
		end

		if root2 < x0 || root2 > x1
			val3 = typemax(T)
		end
	end

	theMin = min(val1, val2, val3, val4)

	minPos = typemax(T)
	if val1 == theMin
		minPos = x0
	elseif val2 == theMin
		minPos = root1
	elseif val3 == theMin
		minPos = root2
	elseif val4 == theMin
		minPos = x1
	else
		error("Failed to find the minimum. This is probably a bug in Splines.jl")
	end

	return (theMin, minPos)
end

function mincubic{T}(x::T, a::T, b::T, c::T, d::T) #search between 0 and x
	return mincubic(zero(T), x, a, b, c, d)
end

function integrate_cubic{T}(x0::T, x1::T, a::T, b::T, c::T, d::T) #a+bx+cx^2+dx^3
  return @evalpoly(x1, 0, a, b/2, c/3, d/4)-@evalpoly(x0, 0, a, b/2, c/3, d/4)
end
