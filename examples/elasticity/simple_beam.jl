using SimpleFiniteElements
import SimpleFiniteElements.Elasticity: ∫∫f_dot_v!, ∫∫λ_div_u_div_v!, 
					∫∫2μ_εu_εv!, elasticity_soln, 
                                        visualise_soln
import StaticArrays: SA
using PyPlot

path = joinpath("..", "..", "spatial_domains", "beam1.geo")
gmodel = GeometryModel(path)
λ = 10.0
μ = 0.5
bilinear_forms = Dict("Beam" => [(∫∫λ_div_u_div_v!, λ),
		                 (∫∫2μ_εu_εv!, μ)])

f = SA[0.0, -0.5]
linear_funcs = Dict("Beam" => (∫∫f_dot_v!, f))
gD = SA[0.0, 0.0]
essential_bcs = [("Fixed", gD)]

hmax = 0.2
mesh = FEMesh(gmodel, hmax)
dof = DegreesOfFreedom(mesh, essential_bcs)

u1h, u2h = elasticity_soln(dof, bilinear_forms, linear_funcs)

scale = 150.0
visualise_soln(dof, u1h, u2h, scale)


