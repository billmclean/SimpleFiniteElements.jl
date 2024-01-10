# Solve Steklov eigenproblem without forming the Schur complement.
using SimpleFiniteElements
import SimpleFiniteElements.Poisson: ∫∫a_∇u_dot_∇v!, ∫c_u_v!
import SimpleFiniteElements.Utils: gmsh2pyplot
import IterativeSolvers: lobpcg
import Printf: @sprintf
using PyPlot

surface = true # true for surface plots, false for contour plots
nev = 3
path = joinpath("..", "..", "spatial_domains", "roundL.geo")
gmodel = GeometryModel(path)
essential_bcs = [ ("Left", 0.0), ("Bottom", 0.0) ]
LHS_bilinear_forms = Dict("Omega" => (∫∫a_∇u_dot_∇v!, 1.0))
RHS_bilinear_forms = Dict("Top"   => (∫c_u_v!, 1.0),
 	                  "Right" => (∫c_u_v!, 1.0))

hmax = 0.05
mesh = FEMesh(gmodel, hmax)
dof = DegreesOfFreedom(mesh, essential_bcs)

A_free, A_fix = assemble_matrix(dof, LHS_bilinear_forms)
B_free, B_fix = assemble_matrix(dof, RHS_bilinear_forms)

largest = false
results = lobpcg(A_free, B_free, largest, nev; maxiter=400)
println("Residual norms:")
display(results.residual_norms)
ϕ_free = results.X
ϕ_fix = zeros(dof.num_fixed)

x, y, triangles = gmsh2pyplot(dof)
for k = 1:nev
    ϕh = [ ϕ_free[:,k]; ϕ_fix ]
    figure(k)
    if surface
        plot_trisurf(x, y, triangles, ϕh; cmap="hsv")
	zlabel("ϕₕ")
    else
        tricontourf(x, y, triangles, ϕh)
        colorbar()
    end
    xlabel("x")
    ylabel("y")
    eigval = @sprintf("%0.4f", results.λ[k])
    s = latexstring("\$\\lambda_", k, "=", eigval, "\$")
    title(s)
end
