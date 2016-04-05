#Test for integral
seed = rand(Int)
print_with_color(:green, "Testing integration")
print_with_color(:blue, "    using seed $seed\n")

max_len = 100
len = 2+abs(rand(Int))%(max_len-2)
knots = Float64[]

while length(knots) < len
	knots = sort(unique([(rand()*len)::Float64 for i in 1:len]))
end

values = [rand() for i in knots]

s = Spline(knots, values)

@test_approx_eq integral(s, knots[1], knots[end]) -1*integral(s, knots[end], knots[1])

k0 = rand(Int)%(len-1)
k1 = k0+1
x0 = s.knots[k0]
x1 = s.knots[k1]

@test_approx_eq integral(s, x0, x1) Splines.integrate_cubic(0.0, x1-x0, s.a[k0], s.b[k0], s.c[k0], s.d[k0])
