
#TODO: write tests for mincubic and maxcubic

function randfloat()
	return 2*(1-rand(Float64))::Float64
end

max_len = 100.0
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
knotsNaN[Int(floor(rand()*len))] = NaN #FIXME: this could generate errors with accessing element 0
valuesNAN[Int(floor(rand()*len))] = NaN #FIXME: this could generate errors with accessing element 0
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
@test domain(s) == (domainMin, domainMax)

value_indomain = rand(domainMin:domainMax)
value_notindomain = value_indomain + domainMax + rand(Float64)
@assert value_notindomain < domainMin || value_notindomain > domainMax

@test_throws DomainError s(value_notindomain)

for (i, knot) in enumerate(knots)
	@test_approx_eq_eps s(knot) values[i] values[i]/10000
end

function test_values()
	tol = rand(Float64)/10
	disc = uniform_discretize(s, tol)
	mesh = Splines.uniform_discretize_mesh(s, tol)
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
#inserting a value actually adds a knot
newknot = Float64(randindex()+randfloat())
newvalue = rand()
insert!(s, newknot, newvalue)
@test s(newknot) == newvalue

#the array version of insert also works
oldknots = copy(s.knots)
oldvalues = copy(s.values)
newknots = [Float64(rand(domain(s)[1]:rand():domain(s)[2])) for i in 1:rand(1:10)]
newvalues = [rand() for i in newknots]
insert!(s, newknots, newvalues)
for (i, knot) in enumerate(newknots)
	@test_approx_eq s(knot) newvalues[i]
end

#test that the new knots were properly added
@test s.knots  == unique(sort(union(oldknots, newknots)))

#test that the previous values are unchanged
for (i, knot) in enumerate(oldknots)
	@test_approx_eq s(knot) oldvalues[i]
end

##### End insert! tests
