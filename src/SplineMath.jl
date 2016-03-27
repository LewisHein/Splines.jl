import Base.Math: +, -, *, /
import Base: abs2

#generic case to apply a binary operation to two splines
function apply_binary_op{T}(f::Spline{T}, g::Spline{T}, op::Function)
	fdomain = domain(f)
	gdomain = domain(g)

	newdomain_low = max(fdomain[1], gdomain[1])
	newdomain_high = min(fdomain[end], gdomain[end])
	newdomain = (newdomain_low, newdomain_high)
	newknots = Array{T, 1}(0)

	for knot in union(f.knots, g.knots)
		if knot >= newdomain_low && knot <= newdomain_high
			push!(newknots, knot)
		end
	end

	newknots = sort(unique(newknots))
	@assert !any(isnan(newknots))

	newvalues = similar(newknots)

	for (i, knot) in enumerate(newknots)
		newvalues[i] = op(f(knot), g(knot))
	end

	return Spline(newknots, newvalues)
end

#Arithmetic for splines is defined just as for functions, since they are really approximation that can live in a finite computer to functions that cannot. Therefore, addition and scalar multiplication are defined just as you would expect

#f(x)+g(x)
function +{T}(f::Spline{T}, g::Spline{T})
	return apply_binary_op(f, g, (x, y)->x+y)
end

function -{T}(f::Spline{T}, g::Spline{T})
	return apply_binary_op(f, g, (x, y)->x-y)
end


#f(x)*g(x)
function *{T}(f::Spline{T}, g::Spline{T})
	return apply_binary_op(f, g, (x, y)->x*y)
end


#f(x)/g(x)
function /{T}(f::Spline{T}, g::Spline{T})
	return apply_binary_op(f,  g, (x, y)->x/y)
end



#f(x)+c
function +{T}(f::Spline{T}, c::Number)
	f_new = deepcopy(f)
	f_new.a += c

	return f_new
end

#c+f(x)
function +{T}(c::Number, f::Spline{T})
	return f+c
end


#c*f(x)
function *{T}(c::Number, f::Spline{T})
	f_new = deepcopy(f)
	f_new.a *= c
	f_new.b *= c
	f_new.c *= c
	f_new.d *= c

	return f_new
end


#f(x)*c
function *{T}(f::Spline{T}, c::Number)
	return c*f
end


#f(x)/c
function /{T}(f::Spline{T}, c::Number)
	f_new = deepcopy(f)
	f_new.a /= c
	f_new.b /= c
	f_new.c /= c
	f_new.d /= c

	return f_new
end

function abs2{T<:Real}(f::Spline{T})
	return f^2
end

###   Now we define function composition for splines and functions, introducing the ∘ operator that does this   ###
function compose{T}(f::Function, g::Spline{T})
	newknots = discretize_mesh(g)
	newvalues = similar(newknots)

	for (i, knot) in enumerate(newknots)
		newvalues[i] = f(g(knot))
	end

	return Spline(newknots, newvalues)
end

function ∘{T}(f::Function, g::Spline{T})
	return compose(f, g)
end

export compose
export ∘
