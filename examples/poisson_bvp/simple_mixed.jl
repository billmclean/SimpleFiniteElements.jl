# Simple mixed BVP.  Use simple_mixed.tex to create simple_mixed.pdf for a 
# description.

using SimpleFiniteElements
import SimpleFiniteElements.Poisson: ∫∫a_∇u_dot_∇v!, ∫∫f_v!, ∫g_v!
import SimpleFiniteElements.Utils: gmsh2pyplot
import LinearAlgebra: cholesky
import IterativeSolvers: cg!
using PyPlot

direct_solver = false
path = joinpath("..", "..", "spatial_domains", "roundL.geo")
gmodel = GeometryModel(path)

bilinear_forms = Dict("Omega" => (∫∫a_∇u_dot_∇v!, 1.0))
linear_functionals = Dict("Omega"  => (∫∫f_v!, 2.0),
                          "Left"   => (∫g_v!, -3/4),
			  "Bottom" => (∫g_v!, -3/4))
essential_bcs = [ ("Top", 0.0), ("Right", 0.0) ]

hmax = 0.15
mesh = FEMesh(gmodel, hmax)
dof = DegreesOfFreedom(mesh, essential_bcs)

A_free, A_fix = assemble_matrix(dof, bilinear_forms)
b_free, u_fix = assemble_vector(dof, linear_functionals)
b = b_free - A_fix * u_fix
if direct_solver
    println("Using sparse Cholesky solver.")
    C = cholesky(A_free)
    u_free = C \ b
else
    println("Using conjugate gradients.")
    u_free = zeros(dof.num_free)
    u_free, ch = cg!(u_free, A_free, b, log=true)
end
uh = [ u_free; u_fix ]

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

if !direct_solver
    residuals = ch.data[:resnorm]
    figure(4)
    semilogy(range(1,ch.iters, length=ch.iters), residuals, "o")
    grid(true)
    xlabel("iterations")
    ylabel("residual norm")
    title("CG iterations")
end
