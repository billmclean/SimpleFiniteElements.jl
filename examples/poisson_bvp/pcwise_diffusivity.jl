# Poisson problem with piecewise-constant diffusivity

using SimpleFiniteElements
import SimpleFiniteElements.Poisson: ∫∫a_∇u_dot_∇v!, ∫∫f_v!, poisson_soln
import SimpleFiniteElements.Utils: gmsh2pyplot
using PyPlot

path = joinpath("..", "..", "spatial_domains", "circle_in_square.geo")
gmodel = GeometryModel(path)
bilinear_forms = Dict("Inner Domain" => (∫∫a_∇u_dot_∇v!, 3.0), 
                      "Outer Domain" => (∫∫a_∇u_dot_∇v!, 1.0))
linear_funcs = Dict("Inner Domain" => (∫∫f_v!, (x,y)->2+x-y),
                    "Outer Domain" => (∫∫f_v!, 0.0))
essential_bcs = [ ("Outer Bdry", 0.0) ]

hmax = 0.25
mesh = FEMesh(gmodel, hmax)
dof = DegreesOfFreedom(mesh, essential_bcs)
uh = poisson_soln(dof, bilinear_forms, linear_funcs)

x, y, triangles = gmsh2pyplot(dof)

figure(1)
triplot(x, y, triangles)
axis("equal")

figure(2)
tricontourf(x, y, triangles, uh, 15)
grid(true)
colorbar()
axis("equal")

figure(3)
plot_trisurf(x, y, triangles, uh, cmap="autumn")
xlabel("x")
ylabel("y")
zlabel("u")

