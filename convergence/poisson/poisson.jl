#
# Simple example of a Poisson problem:
#
#           y
#             |
#             |           u_y = 0
#          Ly ----------------------------------
#             |                                |
#             |                                |
#             |                                |
#       u = 0 |    u - ∇²lxu_xx - u_yy = f         | u = 0
#             |                                |
#             |                                |
#             |                                |
#           -----------------------------------|-----  x
#                         u_y = 0               Lx
#
# Here, 
#
#      u = sin(omegax x) * cos(omegay y) 
#      f = ( 1 + omegax^2 + omegay^2 ) * sin(omegax x) * cos(omegay y)
#      omegax = nx pi / Lx
#      omegay = ny pi / Ly
#

using SimpleFiniteElements
import SimpleFiniteElements.Poisson: ∫∫a_∇u_dot_∇v!, ∫∫c_u_v!, ∫∫f_v!
import SimpleFiniteElements.NonConformingPoisson as NCP
using Printf

include("params.jl")

const ωx = 3π / Lx
const ωy =  π / Ly

function exact_u(x, y)
    return sin(ωx*x) * cos(ωy*y) 
end

function f(x, y)
    c = 1 + ωx^2 + ωy^2
    return c * exact_u(x, y)
end

conforming = false
path = joinpath("..", "..", "spatial_domains", "rect.geo")
gmodel = GeometryModel(path)
println("Poisson problem with zero mixed boundary data.")
if conforming
    println("Using conforming elements.")
    bilinear_forms = [("Omega", ∫∫a_∇u_dot_∇v!, 1.0),
		      ("Omega", ∫∫c_u_v!, 1.0)]
    linear_funcs = [("Omega", ∫∫f_v!, f)]
else
    println("Using non-conforming elements.")
    bilinear_forms = [("Omega", NCP.∫∫a_∇u_dot_∇v!, 1.0),
		      ("Omega", NCP.∫∫c_u_v!, 1.0)]
    linear_funcs = [("Omega", NCP.∫∫f_v!, f)]
end
essential_bcs = [("Left", 0.0), ("Right", 0.0)]

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

