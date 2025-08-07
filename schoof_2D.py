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
v = df.Function(V)  # initialized to zero by default
f = df.Function(V).interpolate(f0 + x/2)
tau_c = df.conditional(abs(x) > 1/2, df.conditional(abs(x) < 3/2, tau0, 0.0), 0.0)  # UFL expression

# 3D
# s = df.sqrt( u.dx(0)**2 +  u.dx(0)*v.dx(1) + v.dx(1)**2 \
#             + 0.25*(u.dx(1)+v.dx(0))**2 + 0.25*u.dx(2)**2 + 0.25*v.dx(2)**2 )

s  = df.sqrt(v.dx(0)**2 + 0.25*v.dx(1)**2 + 1e-15)
G  = 2*B*s**p / p

# Coulomb; id=1 is bed of glacier
#j = tau_c * abs(v) * df.ds(1)    # unregularized
j = tau_c * df.sqrt(v*v + 1.0e-8) * df.ds(1)    # regularized

# objective
J = (G - f * v ) * df.dx + j
F = df.derivative(J, v)

# solver parameters to give some feedback
sp = {
    "snes_linesearch_type": "bt",  # I'm confused about what exactly is needed for linesearch
    "snes_max_it": 200,
    #"snes_rtol": 1.0e-8,
    #"snes_atol": 1.0e-12,
    #"snes_stol": 0.0,
    #"ksp_type": "preonly",
    #"pc_type": "lu",
    #"pc_factor_mat_solver_type": "mumps",
    "snes_monitor": None,
    "snes_converged_reason": None,
}

# solve; note that options like """-s_snes_linesearch_type basic""" can be used
df.solve(F == 0, v, solver_parameters=sp, options_prefix="s")  # add back Dirichlet BC on small part of domain boundary?
v.rename("horizontal velocity")
VTKFile('result_schoof.pvd').write(v)
