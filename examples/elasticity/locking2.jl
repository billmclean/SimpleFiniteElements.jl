using SimpleFiniteElements
import SimpleFiniteElements.Elasticity: ∫g_dot_v!, ∫∫λ_div_u_div_v!, 
                                        ∫∫2μ_εu_εv!, elasticity_soln, visualise_soln
import Printf: @printf
import StaticArrays: SA
import SimpleFiniteElements.Utils: gmsh2pyplot
using PyPlot

path = joinpath("..", "..", "spatial_domains", "cooks_membrane.geo")
gmodel = GeometryModel(path)

E = 1.0
#ν = 0.49999
ν = 0.49
λ = E * ν / ( (1 + ν) * (1 - 2ν) )
μ = E / ( 2 * (1 + ν) )
@printf("Cooks membrane problem with E = %g, λ = %g, μ = %g.\n", E, λ, μ)
gD = SA[0.0, 0.0]
essential_bcs = [ ("Left", gD) ]
c = 0.1
gN = SA[0.0, c]
hmax = 4.0
bilinear_forms = Dict("Omega" => [(∫∫λ_div_u_div_v!, λ),
                                  (∫∫2μ_εu_εv!, μ)])
linear_funcs = Dict("Right" => (∫g_dot_v!, gN))
mesh = FEMesh(gmodel, hmax; order=1, save_msh_file=false,
              refinements=0, verbosity=2)

dof = DegreesOfFreedom(mesh, essential_bcs)
@printf("Number of degrees of freedom = %d\n", 2dof.num_free)
u1h, u2h = elasticity_soln(dof, bilinear_forms, linear_funcs)

scale = 5.0
visualise_soln(dof, u1h, u2h, scale)

x, y, triangles = gmsh2pyplot(dof)
figure(3)
triplot(x, y, triangles)
