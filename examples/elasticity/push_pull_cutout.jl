using SimpleFiniteElements
import SimpleFiniteElements.Elasticity: ∫g_dot_v!, ∫∫λ_div_u_div_v!, 
					∫∫2μ_εu_εv!, elasticity_soln, 
                                        visualise_soln
import StaticArrays: SA
using PyPlot

push = false
path = joinpath("..", "..", "spatial_domains", "cutout2.geo")
gmodel = GeometryModel(path)
λ = 10.0
μ = 0.5
bilinear_forms = Dict("Omega" => [(∫∫λ_div_u_div_v!, λ),
		                  (∫∫2μ_εu_εv!, μ)])

if push
    push_gN(x,y) = SA[-1.0, 0.0]
    linear_funcs = Dict("Traction" => (∫g_dot_v!, push_gN))
else
    pull_gN = SA[1.0, 0.0]
    linear_funcs = Dict("Traction" => (∫g_dot_v!, pull_gN))
end
gD = SA[0.0, 0.0]
essential_bcs = [("Fixed", gD)]

hmax = 0.25
mesh = FEMesh(gmodel, hmax)
dof = DegreesOfFreedom(mesh, essential_bcs)

u1h, u2h = elasticity_soln(dof, bilinear_forms, linear_funcs)

scale = 10.0
visualise_soln(dof, u1h, u2h, scale)
