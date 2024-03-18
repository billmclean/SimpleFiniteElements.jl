using SimpleFiniteElements
import StaticArrays: SA
import SimpleFiniteElements.Elasticity: ∫∫λ_div_u_div_v!, ∫∫2μ_εu_εv!, 
					∫∫f_dot_v!, elasticity_soln, 
                                        L2error 
import SimpleFiniteElements.NonConformingElasticity as NCE
import Printf: @printf

if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
import Example1
import Example2 
import Example3
import Example4

ex = Example2

nrows = 4
conforming = false
quadrature_level = 1
println(ex.msg)
println("Computing L2 error using $(3 * 4^quadrature_level) "
        * "quadrature points per element.")

let
    if conforming
        println("Using conforming elements.")
	bilinear_forms = Dict("Omega" => [(∫∫λ_div_u_div_v!, ex.λ),
					  (∫∫2μ_εu_εv!, ex.μ)])
        linear_funcs = Dict("Omega" => (∫∫f_dot_v!, ex.f))
        mesh_order = 1
    else
        println("Using non-conforming elements.")
        bilinear_forms = Dict("Omega" => [(NCE.∫∫a_div_u_div_v!, ex.μ_plus_λ),
                                          (NCE.∫∫μ_∇u_colon_∇v!, ex.μ),
					  (NCE.correction!, ex.∇μ)])
	linear_funcs = Dict("Omega" => (NCE.∫∫f_dot_v!, ex.f))
        mesh_order = 2
    end

    function fem_error(mesh::FEMesh)
        if conforming
            dof = DegreesOfFreedom(mesh, ex.essential_bcs)
        else
            dof = DegreesOfFreedom(mesh, ex.essential_bcs, NCE.ELT_DOF)
        end
        u1h, u2h = elasticity_soln(dof, bilinear_forms, linear_funcs)
        u = get_nodal_values(ex.exact_u, dof, 2)

        max_error = 0.0
        for k in eachindex(u1h)
            max_error = max(max_error, hypot(u[k,1] - u1h[k], u[k,2] - u2h[k]))
        end
        if conforming
            L2_error = L2error(u1h, u2h, ex.exact_u, dof, quadrature_level) 
        else
            L2_error = NCE.L2error(u1h, u2h, ex.exact_u, dof, quadrature_level) 
        end
        h = max_elt_diameter(mesh)
        num_free = dof.num_free
        return max_error, L2_error, 2num_free, h
    end

    print("Generating $nrows meshes ... ")
    start = time()
    hmax = 0.2
    mesh = FEMesh(ex.gmodel, hmax, order=mesh_order, save_msh_file=false, 
                  refinements=nrows-1, verbosity=2)
    if nrows == 1
	mesh = [mesh]
    end
    elapsed = time() - start
    @printf("in %0.4f seconds.\n\n", elapsed)
    maxerr = Vector{Float64}(undef, nrows)
    L2err = similar(maxerr)
    @printf("%6s  %8s  %12s  %5s  %8s  %5s  %8s\n\n", 
            "h", "unknowns", "nodal error", "rate", "L2 error", "rate", "secs")
    for row in eachindex(maxerr)
        hmax /= 2
        start = time()
        maxerr[row], L2err[row], num_free, h = fem_error(mesh[row])
        elapsed = time() - start
        if row == 1
            @printf("%6.4f  %8d  %12.1e  %5s  %8.1e  %5s  %8.1e\n",
                    h, num_free, maxerr[row], "", L2err[row], "", elapsed)
        else
            maxrate = log2(maxerr[row-1] / maxerr[row])
            L2rate = log2(L2err[row-1] / L2err[row])
 
            @printf("%6.4f  %8d  %12.1e  %5.3f  %8.1e  %5.3f  %8.1e\n", 
                    h, num_free, maxerr[row], maxrate, 
                    L2err[row], L2rate, elapsed)
        end
    end
end
