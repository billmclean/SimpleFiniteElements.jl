# Simple Poisson problem with homogeneous Dirichlet boundary conditions.

using SimpleFiniteElements
import SimpleFiniteElements.Poisson: ∫∫a_∇u_dot_∇v!, ∫∫f_v!
import SimpleFiniteElements.Utils: gmsh2pyplot
import LinearAlgebra: cholesky
using PyPlot

path = joinpath("..", "..", "spatial_domains", "keyhole.geo")
gmodel = GeometryModel(path)
bilinear_forms = Dict("Omega" => (∫∫a_∇u_dot_∇v!, 1.0))
linear_funcs = Dict("Omega" => (∫∫f_v!, 2.0))
essential_bcs = [("Gamma", 2.0)]

hmax = 0.25
mesh = FEMesh(gmodel, hmax)
dof = DegreesOfFreedom(mesh, essential_bcs)

A_free, A_fix = assemble_matrix(dof, bilinear_forms)
b_free, u_fix = assemble_vector(dof, linear_funcs)
b = b_free - A_fix * u_fix
C = cholesky(A_free)
u_free = C \ b
uh = [u_free; u_fix]

x, y, triangles = gmsh2pyplot(dof)

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

