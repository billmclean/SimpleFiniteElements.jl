# Simple example of a Steklov eigenproblem:
#
#           y
#             |
#             |           u_y = λu
#          Ly ----------------------------------
#             |                                |
#             |                                |
#             |                                |
#       u = 0 |        u_xx + u_yy = 0         | u = 0
#             |                                |
#             |                                |
#             |                                |
#           --|--------------------------------|-----  x
#             |            u = 0              Lx
#
# Here, omega = 2 n pi / Lx, and the exact solution is
#
#      u = sin(omega x) * sinh(omega y) / sinh(omega Ly)
#

using SimpleFiniteElements
import SimpleFiniteElements.Poisson: ∫∫a_∇u_dot_∇v!, ∫c_u_v!
import SimpleFiniteElements.NonConformingPoisson as NCP
import LinearAlgebra: eigen
import IterativeSolvers: lobpcg
using Printf

include("params.jl")
nev = 3

const ω = π / Lx

function exact_u(x, y)
    return sin(ω*x) * sinh(ω*y) / sinh(ω*Ly)
end

exact_λ = [ ω * coth(ω*Ly) for ω in [1, 2, 3] .* π / Lx ]
largest = false

conforming = false
path = joinpath("..", "..", "spatial_domains", "rect.geo")
gmodel = GeometryModel(path)
println("Steklov eigenproblem.")
if conforming
    println("Using conforming elements.")
    LHS_bilinear_forms=[("Omega", ∫∫a_∇u_dot_∇v!, 1.0)]
    RHS_bilinear_forms=[("Top", ∫c_u_v!, 1.0)]
else
    println("Using non-conforming elements.")
    LHS_bilinear_forms=[("Omega", NCP.∫∫a_∇u_dot_∇v!, 1.0)]
    RHS_bilinear_forms=[("Top", NCP.∫c_u_v!, 1.0)]
end
essential_bcs = [("Bottom", 0.0), ("Left", 0.0), ("Right", 0.0)]

hmax = h0
err = zeros(nev, refinements+1)
@printf("%5s %10s %8s %10s %8s %10s %8s %8s\n\n", "DoF", 
	"λ₁ error", "rate", "λ₂ error", "rate", "λ₃ error", "rate", "seconds")
for k = 0:refinements
    global hmax
    hmax /= 2
    start = time()
    if conforming
        mesh = FEMesh(gmodel, hmax, order=1, save_msh_file=false)
        dof = DegreesOfFreedom(mesh, essential_bcs)
    else
        mesh = FEMesh(gmodel, hmax, order=2, save_msh_file=false)
        dof = DegreesOfFreedom(mesh, essential_bcs, NCP.ELT_DOF)
    end
    num_steklov = reorder_dof_steklov!(mesh, dof, RHS_bilinear_forms)
    A_free, A_fix = assemble_matrix(mesh, dof, LHS_bilinear_forms)
    B_free, B_fix = assemble_matrix(mesh, dof, RHS_bilinear_forms)
    S = SchurComplement(A_free, num_steklov)
    B11 = B_free[1:num_steklov,1:num_steklov]
    if num_steklov < 40 # use dense eigensolver
	Smat = Matrix(S)
	Bmat = Matrix(B11)
	results = eigen(Smat, Bmat)
	λh = results.values[1:nev]
    else
        S = SchurComplement(A_free, num_steklov)
        results = lobpcg(S, B11, largest, nev; maxiter=400)
        λh = results.λ
    end
    finish = time()
    err[:, k+1] = abs.(λh - exact_λ)
    num_free = dof.num_free
    if k == 0
        @printf("%5d %10.2e %8s %10.2e %8s %10.2e\n", 
		num_free, err[1,k+1], "", err[2,k+1], "", err[3,k+1])
    else
        rate = log2.(err[:,k] ./ err[:,k+1])
        elapsed = finish - start
        @printf("%5d %10.2e %8.4f %10.2e %8.4f %10.2e %8.4f %8.4f\n", 
		num_free, err[1,k+1], rate[1],
		err[2,k+1], rate[2], err[3,k+1], rate[3], elapsed)
    end
end

