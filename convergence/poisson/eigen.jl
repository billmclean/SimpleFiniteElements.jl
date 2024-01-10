# Simple example of an eigenproblem:
#
#           y
#             |
#             |            u_y = 0
#          Ly ----------------------------------
#             |                                |
#             |                                |
#             |                                |
#       u = 0 |     u_xx + u_yy + λ u = 0      | u = 0
#             |                                |
#             |                                |
#             |                                |
#           -----------------------------------|-----  x
#                          u_y = 0              Lx
#
# Here, exact solutions are
#
#      u = sin( n π x / Lx ) cos( k π y / Ly )
#
#                      2               2
#      λ = ( n π / Lx )  + ( k π / Ly )
#
# for n = 1, 2, 3, ... and k = 0, 1, 2, ....

using SimpleFiniteElements
import SimpleFiniteElements.Poisson: ∫∫a_∇u_dot_∇v!, ∫∫c_u_v!
import SimpleFiniteElements.NonConformingPoisson as NCP
using Printf
using IterativeSolvers: lobpcg

include("params.jl")
const nev = 3

λ = Float64[]
for n = 1:nev, k = 0:nev-1
    push!(λ, (n*π/Lx)^2+(k*π/Ly)^2)
end

sort!(λ)

conforming = false
path = joinpath("..", "..", "spatial_domains", "rect.geo")
gmodel = GeometryModel(path)
println("Eigenproblem problem.")
if conforming
    println("Using conforming elements.")
    LHS_bilinear_forms = [("Omega", ∫∫a_∇u_dot_∇v!, 1.0)]
    RHS_bilinear_forms = [("Omega", ∫∫c_u_v!, 1.0)]
    extra_nev = 1
else
    println("Using non-conforming elements.")
    LHS_bilinear_forms = [("Omega", NCP.∫∫a_∇u_dot_∇v!, 1.0)]
    RHS_bilinear_forms = [("Omega", NCP.∫∫c_u_v!, 1.0)]
    extra_nev = 8
end
linear_funcs = Tuple[]
essential_bcs = [("Left", 0.0), ("Right", 0.0)]

hmax = h0
err = zeros(nev, refinements+1)
@printf("Eigenvalue errors\n\n")
@printf("%7s|", "DoF")
for j = 1:nev
    @printf("%14s|", "λ_$j    ")
end 
@printf(" elapsed\n")
@printf("%70s\n", "-"^70)
resid = Float64[]
for k = 0:refinements
    global hmax, resid
    hmax /= 2
    start = time()
    if conforming
        mesh = FEMesh(gmodel, hmax, order=1, save_msh_file=false)
        dof = DegreesOfFreedom(mesh, essential_bcs)
    else
        mesh = FEMesh(gmodel, hmax, order=2, save_msh_file=false)
        dof = DegreesOfFreedom(mesh, essential_bcs, NCP.ELT_DOF)
    end
    A_free, A_fix = assemble_matrix(mesh, dof, LHS_bilinear_forms)
    B_free, B_fix = assemble_matrix(mesh, dof, RHS_bilinear_forms)
    largest = false
    results = lobpcg(A_free, B_free, largest, nev+extra_nev)
    λh = results.λ[1:nev]
    finish = time()
    err[:,k+1] = abs.(λh-λ[1:nev])
    num_free = dof.num_free
    if k == 0
        elapsed = finish - start
        @printf("%7d|", num_free) 
        for j = 1:nev
            @printf("%9.2e %4s|", err[j,k+1], "")
        end
        @printf(" %7.4f\n", elapsed)
    else
        elapsed = finish - start
        @printf("%7d|", num_free)
        for j = 1:nev
            rate = log2(err[j,k]/err[j,k+1])
            @printf("%9.2e %4.2f|", err[j,k+1], rate)
        end
        @printf(" %7.4f\n", elapsed)
    end
end

