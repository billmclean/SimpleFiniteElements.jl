# Simple BVP using nonconforming elements  

using SimpleFiniteElements
import SimpleFiniteElements.NonConformingPoisson as nc
import SimpleFiniteElements.FEM: assemble_matrix, assemble_vector
import LinearAlgebra: cholesky
using PyPlot

path = joinpath("..", "..", "spatial_domains", "keyhole.geo")
gmodel = GeometryModel(path)
bilinear_forms = Dict("Omega" => (nc.∫∫a_∇u_dot_∇v!, 1.0))
linear_functionals = Dict("Omega" => (nc.∫∫f_v!, 2.0))
essential_bcs = [("Gamma", 0.0)]

hmax = 0.25
mesh = FEMesh(gmodel, hmax, order=2)
dof = DegreesOfFreedom(mesh, essential_bcs, nc.ELT_DOF)

A_free, A_fix = assemble_matrix(dof, bilinear_forms)
b_free, u_fix = assemble_vector(dof, linear_functionals)
b = b_free - A_fix * u_fix
C = cholesky(A_free)
u_free = C \ b
uh_midpt = [u_free; u_fix]

x, y, triangles, uh = nc.gmsh2pyplot(dof, uh_midpt)

figure(1)
plot_trisurf(x, y, triangles, uh, cmap="cool")
xlabel("x")
ylabel("y")
zlabel("u")

