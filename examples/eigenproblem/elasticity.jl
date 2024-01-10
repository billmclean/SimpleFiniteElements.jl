using SimpleFiniteElements
import SimpleFiniteElements.Elasticity: ∫∫λ_div_u_div_v!, ∫∫2μ_εu_εv!, ∫∫c_u_dot_v!,
                          visualise_soln
import IterativeSolvers: lobpcg

Lame_λ = 10.0
Lame_μ = 0.5
nev = 3

path = joinpath("..", "..", "spatial_domains", "cutout2.geo")
gmodel = GeometryModel(path)
LHS_bilinear_forms = Dict("Omega" => [(∫∫λ_div_u_div_v!, Lame_λ),
				      (∫∫2μ_εu_εv!, Lame_μ)])
RHS_bilinear_forms = Dict("Omega" => (∫∫c_u_dot_v!, 5.0))
essential_bcs = [("Fixed", 0.0)]

hmax = 0.25
mesh = FEMesh(gmodel, hmax)
dof = DegreesOfFreedom(mesh, essential_bcs)
A_free, A_fix = assemble_matrix(dof, LHS_bilinear_forms, 2)
B_free, B_fix = assemble_matrix(dof, RHS_bilinear_forms, 2)
ϕ1_fix = ϕ2_fix = zeros(dof.num_fixed)

largest = false
results = lobpcg(A_free, B_free, largest, nev)
λ = results.λ
ϕ_free = results.X
num_free, num_fixed = dof.num_free, dof.num_fixed
k = 2
ϕ1h = [ ϕ_free[1:num_free,k]; ϕ1_fix ]
ϕ2h = [ ϕ_free[num_free+1:2*num_free,k]; ϕ2_fix ]

scale = 2.0
visualise_soln(dof,  ϕ1h,  ϕ2h, scale, 1)
visualise_soln(dof, -ϕ1h, -ϕ2h, scale, 3)
