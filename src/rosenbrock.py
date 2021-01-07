from sympy import symbols, simplify, pprint

x = symbols('x')
y = symbols('y')

rosenbrock = (1 - x)**2 + 10*(y - x**2)**2
print('rosenbrock:')
pprint(rosenbrock)
print('d_rosenbrock_dx:')
pprint(simplify(rosenbrock.diff(x)))
print('d_rosenbrock_dy:')
pprint(simplify(rosenbrock.diff(y)))
print('d2_rosenbrock_dx2:')
pprint(simplify(rosenbrock.diff(x).diff(x)))
print('d2_rosenbrock_dxy:')
pprint(simplify(rosenbrock.diff(x).diff(y)))
print('d2_rosenbrock_dy2:')
pprint(simplify(rosenbrock.diff(y).diff(y)))
