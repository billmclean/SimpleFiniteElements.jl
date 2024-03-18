# Simple example of a Neumann problem:
#
#           y
#             |
#             |        ∂u/∂y = cos(ωx)
#         L_y |---------------------------------
#             |                                |
#             |                                |
#             |                                |
#   ∂u/∂x = 0 |         -∇²u + u = 0           | ∂u/∂x = 0
#             |                                |
#             |                                |
#             |                                |
#           --|--------------------------------|-----  x
#             |           ∂u/∂y = 0           L_x
#
# Here, ω = 2nπ / L_x, and the exact solution is
#
#      u = cos(ωx) cosh(αy) / (α sinh(α L_y)
#
# where α² = 1 + ω².

using SimpleFiniteElements
import SimpleFiniteElements.Poisson: ∫∫a_∇u_dot_∇v!, ∫∫c_u_v!, ∫g_v!
import SimpleFiniteElements.NonConformingPoisson as NCP
using Printf

include("params.jl")

const ω = 2π / Lx
const α = sqrt(1 + ω^2)

function exact_u(x, y)
    return cos(ω*x) * cosh(α*y) / (α*sinh(α*Ly))
end

function g(x, y)
    return cos(ω*x)
end

conforming = false
path = joinpath("..", "..", "spatial_domains", "rect.geo")
gmodel = GeometryModel(path)
println("Pure Neumann problem.")
if conforming
    println("Using conforming elements.")
    linear_funcs = [("Top", ∫g_v!, g)]
    bilinear_forms = Dict("Omega" => [(∫∫a_∇u_dot_∇v!, 1.0),
				      (∫∫c_u_v!, 1.0)])
    linear_funcs = Dict("Top" => (∫g_v!, g))
else
    println("Using non-conforming elements.")
    bilinear_forms = Dict("Omega" => [(NCP.∫∫a_∇u_dot_∇v!, 1.0),
				      (NCP.∫∫c_u_v!, 1.0)])
    linear_funcs = Dict("Top" => (NCP.∫g_v!, g))
end
essential_bcs = Tuple[]

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

