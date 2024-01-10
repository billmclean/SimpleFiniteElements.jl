# An eigenproblem described in the LaTeX file eigen.tex,
# using a nonconforming method.

using SimpleFiniteElements
import SimpleFiniteElements.NonConformingPoisson as nc
using IterativeSolvers: lobpcg
import Printf: @sprintf
using PyPlot

surface = true # true for Surface plots; false for contour plots.
nev = 3         # number of eigenvalues

path = joinpath("..", "..", "spatial_domains", "roundL.geo")
gmodel = GeometryModel(path)
LHS_bilinear_forms = Dict("Omega" => (nc.∫∫a_∇u_dot_∇v!, 1.0))
RHS_bilinear_forms = Dict("Omega" => ( nc.∫∫c_u_v!, 1.0))
essential_bcs = [("Top", 0.0), ("Right", 0.0)]

hmax = 0.25
mesh = FEMesh(gmodel, hmax, order=2)
dof = DegreesOfFreedom(mesh, essential_bcs, nc.ELT_DOF)

A_free, A_fix = assemble_matrix(dof, LHS_bilinear_forms)
B_free, B_fix = assemble_matrix(dof, RHS_bilinear_forms)
ϕ_fix = zeros(dof.num_fixed)

largest = false
results = lobpcg(A_free, B_free, largest, nev)
λ = results.λ
ϕ_free = results.X
for k = 1:nev
    local x, y, triangles
    ϕh_midpt = [ ϕ_free[:,k]; ϕ_fix ]
    x, y, triangles, ϕh = nc.gmsh2pyplot(dof, ϕh_midpt)
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
    eigval = @sprintf("%0.4f", λ[k])
    s = latexstring("\$\\lambda_", k, "=", eigval, "\$")
    title(s)
end

