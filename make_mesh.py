import gmsh
import numpy as np

def h(x):
    return 1- x**2/4

x = np.arange(-2,2.1,0.1)
z_surf = h(x)
z_bed  = x[1:-1] * 0

xpts = np.concatenate([x, np.flip(x[1:-1])])
ypts = np.concatenate([z_surf, np.flip(z_bed)])

# generate mesh with gmsh
gmsh.initialize()
geometry = gmsh.model.geo

lc = 0.1
points = [geometry.add_point(xi,yi,0,lc) for (xi,yi) in zip(xpts,ypts)]
lines  = [geometry.add_line(pt1, pt2) for (pt1,pt2) in zip(points, np.concatenate([points[1:],[points[0]]])) ]

face  = geometry.add_curve_loop(lines)
plane = geometry.add_plane_surface([face])

i = len(x) - 1
physical_line = geometry.add_physical_group(1, lines[-i:]) # bed
physical_line = geometry.add_physical_group(1, lines[0:i]) # surface
physical_surface = geometry.add_physical_group(2, [plane])

geometry.synchronize()
gmsh.model.mesh.generate(2)
gmsh.write("schoof.msh")
gmsh.finalize()


# plot
import firedrake as df
from firedrake.pyplot import tripcolor, triplot
import matplotlib.pyplot as plt
mesh = df.Mesh("schoof.msh")
fig, axes = plt.subplots()
triplot(mesh, axes=axes)
axes.legend()
axes.axis("equal")
plt.savefig("schoof_mesh.jpg")
