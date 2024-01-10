# Simple mixed BVP.  Use simple_mixed.tex to create simple_mixed.pdf for a 
# description.

using SimpleFiniteElements
import SimpleFiniteElements.NonConformingPoisson as nc
import LinearAlgebra: cholesky
using PyPlot

path = joinpath("..", "..", "spatial_domains", "roundL.geo")
gmodel = GeometryModel(path)
bilinear_forms = Dict("Omega" => (nc.∫∫a_∇u_dot_∇v!, 1.0))
linear_functionals = Dict("Omega"  => (nc.∫∫f_v!, 2.0),
                          "Left"   => (nc.∫g_v!, -3/4),
		          "Bottom" => (nc.∫g_v!, -3/4))
essential_bcs = [ ("Top", 0.0), ("Right", 0.0) ]

hmax = 0.15
mesh = FEMesh(gmodel, hmax, order=2)
dof = DegreesOfFreedom(mesh, essential_bcs, nc.ELT_DOF)

A_free, A_fix = assemble_matrix(dof, bilinear_forms)
b_free, u_fix = assemble_vector(dof, linear_functionals)
b = b_free - A_fix * u_fix
C = cholesky(A_free)
u_free = C \ b
uh_midpt = [ u_free; u_fix ]

x, y, triangles, uh = nc.gmsh2pyplot(dof, uh_midpt)

figure(1)
triplot(x, y, triangles)
axis("equal")

figure(2)
tricontourf(x, y, triangles, uh)
grid(true)
colorbar()
axis("equal")

figure(3)
plot_trisurf(x, y, triangles, uh, cmap="cool")
xlabel("x")
ylabel("y")
zlabel("u")

