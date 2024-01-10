using SimpleFiniteElements
import SimpleFiniteElements.NonConformingElasticity: ∫∫f_dot_v!, 
    ∫∫a_div_u_div_v!, ∫∫μ_∇u_colon_∇v!, correction!, ELT_DOF,
    elasticity_soln, visualise_soln
import StaticArrays: SA
using PyPlot

path = joinpath("..", "..", "spatial_domains", "rect.geo")
gmodel = GeometryModel(path)
λ = 10.0
μ = 0.5
bilinear_forms = Dict("Omega" => [(∫∫a_div_u_div_v!, μ+λ),
                                 (∫∫μ_∇u_colon_∇v!, μ)])

f = SA[0.0, -15.0]
linear_funcs = Dict("Omega" =>  (∫∫f_dot_v!, f))
gD = SA[0.0, 0.0]
essential_bcs = [("Bottom", gD), ("Right", gD), ("Top", gD), ("Left", gD)]

hmax = 0.1
mesh = FEMesh(gmodel, hmax; order=2)
dof = DegreesOfFreedom(mesh, essential_bcs, ELT_DOF)

u1h, u2h = elasticity_soln(dof, bilinear_forms, linear_funcs)

scale = 1.0
edenum = EdgeEnumeration(mesh)
visualise_soln(dof, edenum, u1h, u2h, scale)


