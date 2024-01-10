using SimpleFiniteElements
import SimpleFiniteElements.Elasticity: fundamental_soln, 
           linear_elasticity_soln, visualise_soln, get_nodal_values,
           ∫∫λ_div_u_div_v!, ∫∫2μ_εu_εv!
import SimpleFiniteElements.Utils: gmsh2pyplot
using StaticArrays: SA
using PyPlot

gmodel = GeometryModel("../examples/spatial_domains/cutout1.geo")
λ = 1.0
μ = 0.5
function g(x, y) 
    G = fundamental_soln(SA[x, y], SA[0.0, 1.2], λ, μ)
    return G * SA[0.0, -10.0]
end
bilinear_forms = [("Omega", ∫∫λ_div_u_div_v!, λ),
                  ("Omega", ∫∫2μ_εu_εv!, μ)]
linear_funcs = Tuple[]
essential_bcs = [("Gamma", g)]

hmax = 0.2
mesh = FEMesh(gmodel, hmax)
dof = DegreesOfFreedom(mesh, essential_bcs)

b_free, u_fix = assemble_vector(mesh, dof, linear_funcs, 2)
A_free, A_fix = assemble_matrix(mesh, dof, bilinear_forms, 2)
b = b_free - A_fix * u_fix
u_free = A_free \ b
num_free, num_fixed = dof.num_free, dof.num_fixed
u1h = [ u_free[1:num_free]; u_fix[1:num_fixed] ]
u2h = [ u_free[num_free+1:2*num_free]; u_fix[num_fixed+1:2*num_fixed] ]

#u1h, u2h = linear_elasticity_soln(mesh, dof, bilinear_forms, linear_funcs, 
#				  essential_bcs)

scale = 6.0
visualise_soln(mesh, dof, u1h, u2h, scale)

x, y, triangles = gmsh2pyplot(mesh, dof)
u1, u2 = get_nodal_values(g, mesh, dof)

figure(3)
tricontourf(x, y, triangles, u1h)
axis("equal")
colorbar()
title(L"Finite element approximation to $u_1$")

figure(4)
tricontourf(x, y, triangles, u1)
axis("equal")
colorbar()
title(L"The analytical solution $u_1$")

figure(5)
tricontourf(x, y, triangles, u2h)
axis("equal")
colorbar()
title(L"Finite element approximation to $u_2$")

figure(6)
tricontourf(x, y, triangles, u2)
axis("equal")
colorbar()
title(L"The analytical solution $u_2$")
