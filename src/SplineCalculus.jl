function derivative{T}(s::Spline{T})
	nSplines = length(s.a)

	sDeriv = deepcopy(s)
	for i in 1:nSplines
		sDeriv.values[i] = sDeriv.b[i]
		sDeriv.a[i] = sDeriv.b[i]
		sDeriv.b[i] = sDeriv.c[i]*2
		sDeriv.c[i] = sDeriv.d[i]*3
		sDeriv.d[i] = 0
	end
	sDeriv.values[nSplines+1] = @evalpoly(s.knots[nSplines+1]-s.knots[nSplines], sDeriv.a[end], sDeriv.b[end], sDeriv.c[end], sDeriv.d[end]) #TODO: not tested

	return sDeriv
end

function derivative2{T}(s::Spline{T})
	nSplines = length(s.a)
	sDeriv = deepcopy(s)
	for i in 1:nSplines
		sDeriv.a[i] = sDeriv.c[i]*2
		sDeriv.b[i] = sDeriv.d[i]*6
		sDeriv.c[i] = zero(T)
		sDeriv.d[i] = zero(T)
	end
	sDeriv.values[nSplines+1] = 0
	return sDeriv
end

"""Compute the maximum absolute value of the second derivative of s"""
function maxabs2deriv(s::Spline)
	maxA2D = sum = 0
	for i in 1:length(s.knots)-1
		maxA2D = max(maxA2D, abs(sum += (2*s.c[i]+6*s.d[i]*(s.knots[i+1]-s.knots[i]))))
	end

	return maxA2D
end

function integral{T}(s::Spline{T}, a::T, b::T)
  if !indomain(s, a) || !indomain(s, b)
    throw(DomainError())
  end
  if b < a
    return -1*integral(s, b, a)
  end
  if b == a
    return 0
  end

  startknot = findlast(x->x<=a, s.knots)
  endknot = findlast(x->x<b, s.knots)

	if startknot == endknot
		integral_val = integrate_cubic(a-s.knots[startknot], b-s.knots[startknot], s.a[startknot], s.b[startknot], s.c[startknot], s.d[startknot])
	else
  	integral_val = zero(T)
  	for knot in startknot:endknot
    	integral_val += integrate_cubic(0.0, s.knots[knot+1]-s.knots[knot], s.a[knot], s.b[knot], s.c[knot], s.d[knot])
  	end
  	integral_val -= integrate_cubic(0.0, a-s.knots[startknot], s.a[startknot], s.b[startknot], s.c[startknot], s.d[startknot])
  	integral_val -= integrate_cubic(b-s.knots[endknot], s.knots[endknot+1]-s.knots[endknot], s.a[endknot], s.b[endknot], s.c[endknot], s.d[endknot])
	end

  return integral_val
end

function integral_squared{T}(s::Spline{T}, a::T, b::T)
	return integral(s^2, a, b)
end

export derivative
export derivative2
export maxabs2deriv
export integral
export integral_squared
