#
# Simple example of a Poisson problem with a variable coefficient:
#
#           y
#             |
#             |           u_y = 0
#          Ly ----------------------------------
#             |                                |
#             |                                |
#             |                                |
#       u = 0 |   -(au_x)_x - (au_y)_y = f     | u = 0
#             |                                |
#             |                                |
#             |                                |
#           -----------------------------------|-----  x
#                         u_y = 0               Lx
#
# Here, 
#
#      u = sin(ωx x) * cos(ωy y) 
#      a = exp(x-y)
#      f = ω exp(x-y) ( 2ω sin(ω x) cos(ω y) - cos ω(x-y) )
#      ωx = mπ / Lx
#      ωy = nπ / Ly
#

using SimpleFiniteElements
import SimpleFiniteElements.Poisson: ∫∫a_∇u_dot_∇v!, ∫∫f_v!, error_norms
import SimpleFiniteElements.NonConformingPoisson as NCP
import StaticArrays: SA
using Printf

include("params.jl")

const ωx = π / Lx
const ωy = π / Ly

exact_u(x, y) = sin(ωx*x) * cos(ωy*y)
∇u(x, y) = SA[ ωx * cos(ωx*x) * cos(ωy*y), 
	      -ωy * sin(ωx*x) * sin(ωy*y)]

a(x, y) = exp(x-y)

f(x, y) = exp(x-y) * ( 
             ωx * ( ωx * sin(ωx*x) - cos(ωx*x) ) * cos(ωy*y)
           + ωy * ( ωy * cos(ωy*y) - sin(ωy*y) ) * sin(ωx*x)
           )

conforming = true
quadrature_level = 0
path = joinpath("..", "..", "spatial_domains", "rect.geo")
gmodel = GeometryModel(path)
println("Poisson problem with variable coefficient.")
println("Geometry from ", path)
println("Computing L2 error using $(3 * 4^quadrature_level) "
        * "quadrature points per element.")
if conforming
    println("Using conforming elements.")
    bilinear_forms = Dict("Omega" => (∫∫a_∇u_dot_∇v!, a))
    linear_funcs = Dict("Omega" => (∫∫f_v!, f))
else
    println("Using non-conforming elements.")
    bilinear_forms = Dict("Omega" => (NCP.∫∫a_∇u_dot_∇v!, a))
    linear_funcs = Dict("Omega" => (NCP.∫∫f_v!, f))
end
essential_bcs = [("Left", 0.0), ("Right", 0.0)]

hmax = h0
maxerr = zeros(refinements+1)
L2err = similar(maxerr)
H1err = similar(maxerr)
@printf("\n%10s  %12s  %8s  %12s  %8s  %8s\n\n", 
        "DoF", "L2 error", "rate", "H1 error", "rate", "seconds")
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
    maxerr[k+1] = maximum(abs.(uh-u))
    if conforming
	L2err[k+1], H1err[k+1] = error_norms(uh, exact_u, ∇u, 
					     dof, quadrature_level)
    else
	L2err[k+1], H1err[k+1] = NCP.error_norms(uh, exact_u, ∇u, 
				  	         dof, quadrature_level)
    end
    finish = time()
    num_free = length(u_free)
    if k == 0
        @printf("%10d  %12.4e  %8s  %12.4e\n", 
                num_free, L2err[k+1], "", H1err[k+1])
    else
        L2rate = log2(L2err[k]/L2err[k+1])
        H1rate = log2(H1err[k]/H1err[k+1])
        elapsed = finish - start
        @printf("%10d  %12.4e  %8.4f  %12.4e  %8.4f  %8.4f\n", 
                num_free, L2err[k+1], L2rate, H1err[k+1], H1rate, elapsed)
    end
end

