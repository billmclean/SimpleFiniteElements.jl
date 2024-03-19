using SimpleFiniteElements
import SimpleFiniteElements.Elasticity: ∫∫f_dot_v!, ∫∫λ_div_u_div_v!, 
					∫∫2μ_εu_εv!, elasticity_soln, 
                                        visualise_soln
import StaticArrays: SA
using PyPlot

path = joinpath("..", "..", "spatial_domains", "rect.geo")
gmodel = GeometryModel(path)
λ = 10.0
μ = 0.5
bilinear_forms = Dict("Omega" => [(∫∫λ_div_u_div_v!, λ),
		                  (∫∫2μ_εu_εv!, μ)])

f = SA[0.0, -15.0]
linear_funcs = Dict("Omega" => (∫∫f_dot_v!, f))
gD = SA[0.0, 0.0]
essential_bcs = [("Bottom", gD), ("Right", gD), ("Top", gD), ("Left", gD)]

hmax = 0.1
mesh = FEMesh(gmodel, hmax)
dof = DegreesOfFreedom(mesh, essential_bcs)

u1h, u2h = elasticity_soln(dof, bilinear_forms, linear_funcs)

scale = 1.0
visualise_soln(dof, u1h, u2h, scale)

