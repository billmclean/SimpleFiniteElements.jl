using SimpleFiniteElements
import SimpleFiniteElements.Elasticity: ∫g_dot_v!, ∫∫λ_div_u_div_v!, 
					∫∫2μ_εu_εv!, elasticity_soln, 
                                        visualise_soln
import StaticArrays: SA
import SimpleFiniteElements.Utils: gmsh2pyplot
using PyPlot

path = joinpath("..", "..", "spatial_domains", "beam2.geo")
gmodel = GeometryModel(path)
λ = 10.0
μ = 1.5
bilinear_forms = [("Beam", ∫∫λ_div_u_div_v!, λ),
		  ("Beam", ∫∫2μ_εu_εv!, μ)]

push_gN = SA[0.0, -0.25]
linear_funcs = [("Traction", ∫g_dot_v!, push_gN)]
gD = SA[0.0, 0.0]
essential_bcs = [("Fixed", gD)]

hmax = 0.25
mesh = FEMesh(gmodel, hmax)
dof = DegreesOfFreedom(mesh, essential_bcs)

u1h, u2h = elasticity_soln(dof, bilinear_forms, linear_funcs)

scale = 50.0
visualise_soln(dof, u1h, u2h, scale)
