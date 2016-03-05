
#TODO: write tests for mincubic and maxcubic

function randfloat()
	return 2*(1-rand(Float64))::Float64
end

max_len = 10000.0
len = rand(1:max_len)
knots = sort(unique([randfloat()*len::Float64 for i in 1:len]))
len = length(knots)

function randindex()
	return rand(1:len)
end

values = [randfloat() for i in knots]

@assert isa(knots, Array{Float64, 1})
@assert isa(values, Array{Float64, 1})


wrongLength = rand(union(1:len-1, len+1:max_len))
values_wrongLength = [randfloat() for i in 1:wrongLength]
@assert length(values_wrongLength) != length(values)

@test_throws ArgumentError Spline(knots, values_wrongLength)
s = Spline(knots, values)

#Test for refusal to work with NaN values
knotsNaN = copy(knots)
valuesNAN = copy(values)
knotsNaN[Int(floor(rand()*len))] = NaN
valuesNAN[Int(floor(rand()*len))] = NaN
@test_throws ArgumentError Spline(knotsNaN, values)
@test_throws ArgumentError Spline(knots, valuesNAN)
@test_throws ArgumentError Spline(knotsNaN, valuesNAN)

#Test for refusal to work with duplicated knots
knotsDup = copy(knots)
dupInd = Int(floor(rand()*len))
dupInd2 = dupInd
while dupInd2 == dupInd
	dupInd2 = Int(floor(rand()*len))
end
knotsDup[dupInd] = knotsDup[dupInd2]
@test_throws ArgumentError Spline(knotsDup, values)


domainMin = minimum(knots)
domainMax = maximum(knots)
@test domain(s) == domainMin:domainMax

value_indomain = rand(domainMin:domainMax)
value_notindomain = value_indomain + domainMax + rand(Float64)
@assert value_notindomain < domainMin || value_notindomain > domainMax

@test_throws DomainError s(value_notindomain)

for (i, knot) in enumerate(knots)
	@test_approx_eq_eps s(knot) values[i] values[i]/10000
end

function test_values()
	tol = rand(Float64)/10
	disc = discretize(s, tol)
	mesh = Splines.discretize_mesh(s, tol)
	for i in eachindex(disc)
		@test s(mesh[i]) ==  disc[i]
	end
end
test_values()

@test_approx_eq s(xmax(s))  maximum(s)

###### Now test the functions that apply to splines #######

##### recalculate!
values2 = [randfloat() for i in knots]
@test_throws ArgumentError Splines.recalculate!(s, valuesNAN)
@test_throws ArgumentError Splines.recalculate!(s, values_wrongLength)
Splines.recalculate!(s, values2)
test_values()
##### End recalculate! tests

##### insert!
# just change value if knot is already there
insert!(s, Float64(randindex()),  rand()) 
test_values()
