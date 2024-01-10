# More complicated mixed BVP.  Use complicated_mixed.tex to create 
# complicated_mixed.pdf for a description.

using SimpleFiniteElements
import SimpleFiniteElements.Poisson: ∫∫a_∇u_dot_∇v!, ∫∫f_v!, ∫g_v!
import SimpleFiniteElements.Utils: gmsh2pyplot
import LinearAlgebra: factorize
using PyPlot

path = joinpath("..", "..", "spatial_domains", "complicated_keyhole.geo")
gmodel = GeometryModel(path)

bilinear_forms = Dict("LightBlue"  => (∫∫a_∇u_dot_∇v!, 1.0),
		      "LightGreen" => (∫∫a_∇u_dot_∇v!, 10.0))
linear_funcs = Dict("LightBlue"  => (∫∫f_v!, 1.0), 
		    "LightGreen" => (∫∫f_v!, 4.0),
		    "Violet" => (∫g_v!, -1.0))
gD(x, y) = -hypot(x, y) / 2
essential_bcs = [("DarkBlue", gD), ("Black", 0.0)]

hmax = 0.1
mesh = FEMesh(gmodel, hmax)
dof = DegreesOfFreedom(mesh, essential_bcs)

A_free, A_fix = assemble_matrix(dof, bilinear_forms)
b_free, u_fix = assemble_vector(dof, linear_funcs)
b = b_free - A_fix * u_fix
F = factorize(A_free)
u_free = F \ b
uh = [ u_free; u_fix ]

x, y, triangles = gmsh2pyplot(dof)

figure(1)
triplot(x, y, triangles)
axis("equal")

figure(2)
tricontourf(x, y, triangles, uh, 20)
grid(true)
colorbar()
axis("equal")

figure(3)
plot_trisurf(x, y, triangles, uh, cmap="cool")
xlabel("x")
ylabel("y")
zlabel("u")
