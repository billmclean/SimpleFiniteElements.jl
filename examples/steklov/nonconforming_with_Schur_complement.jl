# Solve Steklov eigenproblem via a Schur complement.
using SimpleFiniteElements
import SimpleFiniteElements.MeshGen: reorder_dof_steklov!
import SimpleFiniteElements.NonConformingPoisson as nc
import LinearAlgebra: eigen
import IterativeSolvers: lobpcg
import Printf: @sprintf
using PyPlot

surface = true # true for surface plots, false for contour plots
nev = 3        # number of eigenvalues

path = joinpath("..", "..", "spatial_domains", "roundL.geo")
gmodel = GeometryModel(path)
LHS_bilinear_forms = Dict("Omega" => (nc.∫∫a_∇u_dot_∇v!, 1.0))
RHS_bilinear_forms = Dict("Top"   => (nc.∫c_u_v!, 1.0),
		          "Right" => (nc.∫c_u_v!, 1.0))
essential_bcs = [ ("Left", 0.0), ("Bottom", 0.0) ]

hmax = 0.075
mesh = FEMesh(gmodel, hmax, order=2)
dof = DegreesOfFreedom(mesh, essential_bcs, nc.ELT_DOF)
num_steklov = reorder_dof_steklov!(dof, RHS_bilinear_forms)

A_free, A_fix = assemble_matrix(dof, LHS_bilinear_forms)
B_free, B_fix = assemble_matrix(dof, RHS_bilinear_forms)
B11 = B_free[1:num_steklov,1:num_steklov]

S = SchurComplement(A_free, num_steklov)
if size(S, 1) < 40
    println("Using dense eigensolver")
    Smat = Matrix(S)
    Bmat = Matrix(B11)
    results = eigen(Smat, Bmat)
    ϕ1_free = results.vectors
    λ = results.values
else
    println("Using iterative eigensolver")
    largest = false
    results = lobpcg(S, B11, largest, nev; maxiter=400)
    println("Residual norms:")
    display(results.residual_norms)
    ϕ1_free = results.X
    λ = results.λ
end
A21 = transpose(S.A12)
ϕ2_free = - ( S.CA \ ( A21 * ϕ1_free ) )
ϕ_free = [ ϕ1_free; ϕ2_free ]
ϕ_fix = zeros(dof.num_fixed)

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
