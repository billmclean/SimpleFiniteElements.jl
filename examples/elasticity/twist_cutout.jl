using SimpleFiniteElements
import SimpleFiniteElements.Elasticity: ∫g_dot_v!, ∫∫λ_div_u_div_v!, 
					∫∫2μ_εu_εv!, elasticity_soln, 
                                        visualise_soln
import StaticArrays: SA
using PyPlot

dirichlet = false
path = joinpath("..", "..", "spatial_domains", "cutout2.geo")
gmodel = GeometryModel(path)
λ = 10.0
μ = 0.5
bilinear_forms = Dict("Omega" => [(∫∫λ_div_u_div_v!, λ),
		                  (∫∫2μ_εu_εv!, μ)])

gD = SA[0.0, 0.0]
essential_bcs = Tuple[("Fixed", gD)]
if dirichlet
    gD_twist(x, y) = SA[y, 0.0]
    push!(essential_bcs, ("Traction", gD_twist))
    linear_funcs = Dict()
else
    gN(x, y) = SA[y/3, 0.0]
    linear_funcs = Dict("Traction" => (∫g_dot_v!, gN))
end

hmax = 0.25
mesh = FEMesh(gmodel, hmax)
dof = DegreesOfFreedom(mesh, essential_bcs)

u1h, u2h = elasticity_soln(dof, bilinear_forms, linear_funcs)

if dirichlet
    scale = 5.0
else
    scale = 10.0
end
visualise_soln(dof, u1h, u2h, scale)
