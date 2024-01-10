# Simple example of a Dirichlet problem:
#
#           y
#             |
#             |       u = sin(omega x)
#         L_y |---------------------------------
#             |                                |
#             |                                |
#             |                                |
#       u = 0 |           ∇²u = 0              | u = 0
#             |                                |
#             |                                |
#             |                                |
#           --|--------------------------------|-----  x
#             |            u = 0              L_x
#
# Here, ω = 2nπ / L_x, and the exact solution is
#
#      u = sin(ωx) sinh(ωy) / sinh(ω L_y)
#

using SimpleFiniteElements
import SimpleFiniteElements.Poisson: ∫∫a_∇u_dot_∇v!
import SimpleFiniteElements.NonConformingPoisson as NCP
using Printf

include("params.jl")

const ω = 2π / Lx

exact_u(x, y) = sin(ω*x) * sinh(ω*y) / sinh(ω*Ly)
g(x, y) = sin(ω*x)

conforming = false
path = joinpath("..", "..", "spatial_domains", "rect.geo")
gmodel = GeometryModel(path)
println("Pure Dirichlet problem.")
if conforming
    println("Using conforming elements.")
    bilinear_forms = [("Omega", ∫∫a_∇u_dot_∇v!, 1.0)]
else
    println("Using non-conforming elements.")
    bilinear_forms = [("Omega", NCP.∫∫a_∇u_dot_∇v!, 1.0)]
end
linear_funcs = Tuple[]
essential_bcs = [("Top", g), ("Bottom", 0.0), ("Left", 0.0), ("Right", 0.0)]

hmax = h0
maxerr = zeros(refinements+1)
@printf("%10s  %12s  %8s  %8s\n\n", 
        "DoF", "max error", "rate", "seconds")
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
    A_free, A_fix = assemble_matrix(dof, bilinear_forms)
    b_free, u_fix = assemble_vector(dof, linear_funcs)
    b = b_free - A_fix * u_fix
    u_free = A_free \ b
    uh = [ u_free; u_fix ]
    u = get_nodal_values(exact_u, dof)
    finish = time()
    maxerr[k+1] = maximum(abs.(uh-u))
    num_free = length(u_free)
    if k == 0
        @printf("%10d  %12.4e\n", num_free, maxerr[k+1])
    else
        rate = log2(maxerr[k]/maxerr[k+1])
        elapsed = finish - start
        @printf("%10d  %12.4e  %8.4f  %8.4f\n", 
                num_free, maxerr[k+1], rate, elapsed)
    end
end

