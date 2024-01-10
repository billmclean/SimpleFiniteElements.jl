using SimpleFiniteElements
import SimpleFiniteElements.Elasticity: elasticity_soln, visualise_soln,
                          ∫∫λ_div_u_div_v!, ∫∫μ_εu_εv!, ∫∫f_dot_v!
import SimpleFiniteElements.NonConformingElasticity as nce
import Printf: @printf, @sprintf
import StaticArrays: SA
import PyPlot: figure, title, close

include("Example1.jl")
import .Example1: λ, μ, u, f, msg

α = 0.0
Λ = 1.0
println("α = $α, Λ = $Λ")
scale = 20.0
println(msg)
conforming = false
display_soln = true
gmodel = GeometryModel("../../spatial_domains/square.geo")

if conforming
    bilinear_forms = [("Omega", ∫∫λ_div_u_div_v!, (x, y) ->  λ(x, y, Λ)),
		      ("Omega", ∫∫μ_εu_εv!, (x, y) -> μ(x, y, α))]
    linear_funcs = [("Omega", ∫∫f_dot_v!, (x, y) -> f(x, y, α, Λ))]
    mesh_order = 1
    println("Elasticity with conforming FEM")
else
    bilinear_forms = [("Omega", nce.∫∫λ_div_u_div_v!, (x, y) ->  λ(x, y, Λ)),
		      ("Omega", nce.∫∫μ_εu_εv!, (x, y) -> μ(x, y, α))]
    linear_funcs = [("Omega", nce.∫∫f_dot_v!, (x, y) -> f(x, y, α, Λ))]
    mesh_order = 2
    println("Elasticity with non-conforming FEM")
end

gD = SA[0.0, 0.0]
essential_bcs = [("Top", gD), ("Bottom", gD), ("Left", gD), ("Right", gD)]

if display_soln
    hmax = 0.2
    mesh = FEMesh(gmodel, hmax; order=mesh_order, verbosity=2)
    if conforming
        dof = DegreesOfFreedom(mesh, essential_bcs)
        u1h, u2h = elasticity_soln(mesh, dof, bilinear_forms, linear_funcs)
        u_mat = get_nodal_values(u, mesh, dof, 2)
        visualise_soln(mesh, dof, u1h, u2h, scale)
        visualise_soln(mesh, dof, u_mat[:,1], u_mat[:,2], scale, 3)
        visualise_soln(mesh, dof, u_mat[:,1]-u1h, u_mat[:,2]-u2h, 1.0, 5)
    else
        dof = DegreesOfFreedom(mesh, essential_bcs, nce.ELT_DOF)
        u1h, u2h = elasticity_soln(mesh, dof, bilinear_forms, linear_funcs)
        u_mat = get_nodal_values(u, mesh, dof, 2)
	edenum = EdgeEnumeration(mesh)
        nce.visualise_soln(mesh, dof, edenum, u1h, u2h, scale)
	nce.visualise_soln(mesh, dof, edenum, u_mat[:,1], u_mat[:,2], scale, 3)
    end
    figure(1)
    title("Numerical Solution")
    close(2)
    figure(3)
    title("Exact Solution")
    close(4)
    figure(5)
    title("Error")
    close(6)
end

function fem_error(mesh::FEMesh, verbosity=2)
    if conforming
        dof = DegreesOfFreedom(mesh, essential_bcs)
    else
        dof = DegreesOfFreedom(mesh, essential_bcs, nce.ELT_DOF)
    end
    u1h, u2h = elasticity_soln(mesh, dof, bilinear_forms, linear_funcs)

    u_mat = get_nodal_values(u, mesh, dof, 2)
#    display([u1h[1:10,1] u_mat[1:10,1]])
    err = 0.0
    for k = 1:length(u1h)
        err = max(err, hypot(u_mat[k,1]-u1h[k], u_mat[k,2]-u2h[k]))
    end
    h = max_elt_diameter(mesh)
    return err, 2*dof.num_free, h
end

let
    nrows = 5
    print("Generating ", nrows, " meshes ...")
    start = time()
    hmax = 0.2
    mesh = FEMesh(gmodel, hmax, order=mesh_order, save_msh_file=true,
                  refinements=nrows-1, verbosity=2)
    finish = time()
    elapsed = @sprintf("%0.4f", finish-start)
    println(" in ", elapsed, " seconds.")
    err = Vector{Float64}(undef, nrows)
    @printf("\n%6s  %8s  %8s  %5s  %8s\n\n",
             "h", "unknowns", "error", "rate", "secs")

    for row = 1:nrows
        hmax /= 2
        start = time()
        err[row], num_free, h = fem_error(mesh[row])
        finish = time()
        if row == 1
            @printf("%6.4f  %8d  %8.1e  %5s  %8.1e\n",
                    h, num_free, err[row], "", finish-start)
        else
            rate = log2(err[row-1]/err[row])
            @printf("%6.4f  %8d  %8.1e  %5.3f  %8.1e\n",
                    h, num_free, err[row], rate, finish-start)
        end
    end
end


