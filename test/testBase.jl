
#TODO: write tests for mincubic and maxcubic
seed = rand(Int)
print_with_color(:green, "Testing basic construction")
print_with_color(:blue, "    using seed $seed")
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
knots_somesize = [Float64(i) for i in rand(3:10)]
values_someothersize = [rand() for i in 11:20]
@test_throws ArgumentError insert!(s, knots_somesize, values_someothersize) #test that arrays of different sizes are rejected
oldknots = copy(s.knots)
oldvalues = copy(s.values)
newknots = [Float64(rand(domain(s)[1]:rand():domain(s)[2])) for i in 1:rand(1:10)]

#now ensure that some knots are duplicated to make sure this is handled correctly
for i in 1:rand(2:len)
	insert!(newknots, rand(1:length(newknots)), oldknots[rand(1:len)])
end

newvalues = [rand() for i in newknots]
insert!(s, copy(newknots), copy(newvalues))
for (i, knot) in enumerate(newknots)
#	@test_approx_eq s(knot) newvalues[i]
end

#test that the new knots were properly added
@test s.knots  == unique(sort(union(oldknots, newknots)))

#test that the previous values are unchanged
for (i, knot) in enumerate(oldknots)
	if !(knot in newknots)
		@test_approx_eq s(knot) oldvalues[i]
	end
end

##### End insert! tests

##### moveknot! tests
seed = rand(Int)
print_with_color(:green, "Testing moveknot!")
print_with_color(:blue, "    using seed $seed")
srand(1234)

knots = sort(unique([rand(1.0:1.0:100) for i in 1:10]))
nKnots = length(knots)
values = [rand() for i in knots]
s = Spline(knots, values)
sdomain = domain(s)[1]:rand():domain(s)[2]
@test_throws DomainError moveknot!(s, 1, maximum(knots)+rand()) #can't move outside the domain
@test_throws BoundsError moveknot!(s, -1, rand(sdomain)) #can't move a knot at an out-of-bounds index
@test_throws BoundsError moveknot!(s, nKnots+rand(1:10), rand(sdomain)) #can't move a knot at an out-of-bounds index
@test_throws ArgumentError moveknot!(s, 1, knots[rand(1:nKnots)]) #can't move a knot to where there already is one
knottomove = rand(1:nKnots)
newlocation = knots[knottomove]+rand()
moveknot!(s, knottomove, newlocation)
@test_approx_eq s(newlocation) values[knottomove]

#for the array version of moveknot!
knots = sort(unique([rand(1.0:1.0:100)::Float64 for i in 1:10]))
values = [rand() for i in knots]
s = Spline(knots, values)
newlocations = [(i+rand())::Float64 for i in knots]
newlocations[end] = knots[end]-rand()
moveknot!(s, 1:length(knots), newlocations)
for i in eachindex(newlocations)
	@test_approx_eq s(newlocations[i]) values[i]
end

##### end moveknot! tests
