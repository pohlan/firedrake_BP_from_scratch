import firedrake as df
from firedrake.output import VTKFile

# parameters
p = 4/3
B = 0.1
D0 = 1e-4
tau0 = 0.25   # [0.3, 0.25, 0.2]
f0   = 0      # [0.0, 0.075, 0.15]

# mesh
mesh = df.Mesh('schoof.msh')
x, y = df.SpatialCoordinate(mesh)

V = df.FunctionSpace(mesh, 'CG', 1)
v = df.Function(V)
f = df.Function(V).interpolate(f0 - x/2)
tau_c = df.Function(V).interpolate(df.min_value(df.ge(abs(x),0.5), df.le(abs(x),3/2)) * tau0)

# 3D
# s = df.sqrt( u.dx(0)**2 +  u.dx(0)*v.dx(1) + v.dx(1)**2 \
#             + 0.25*(u.dx(1)+v.dx(0))**2 + 0.25*u.dx(2)**2 + 0.25*v.dx(2)**2 )

s  = df.sqrt(v.dx(0)**2 + 0.25*v.dx(1)**2 + 1e-15)
G  = 2*B*s**p / p

# Coulomb
j = tau_c * abs(v) * df.ds(1)    # id 1 : bed of gllacier

# objective
J = (G - f * v ) * df.dx + j
F = df.derivative(J, v)

# solve
bc = df.DirichletBC(V, 0.0, 'on_boundary')
df.solve(F == 0, v, bcs=[bc,])
VTKFile('result_schoof.pvd').write(v)
