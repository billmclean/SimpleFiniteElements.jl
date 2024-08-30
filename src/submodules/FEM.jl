module FEM

import ..FuncOrConst, ..DegreesOfFreedom, ..GeometryModel, ..FEMesh
import ..MeshGen: elt_node_coord! 
import SparseArrays: sparse, SparseMatrixCSC, nnz
import StaticArrays: SVector

include("deprecated.jl")

function assemble_vector(dof::DegreesOfFreedom,
        linear_functionals::Dict, num_dependent_vars::Int64=1)
    mesh = dof.mesh
    num_free, node_tag_to_dof = dof.num_free, dof.node_tag_to_dof
    # Accumulate the global vector in b_free
    b_free = zeros(num_free * num_dependent_vars)
    for (name, lfuncs) in linear_functionals
	if lfuncs isa Tuple
	    lfuncs = [lfuncs]
	end
	validate_linear_functional(name, lfuncs, dof)
        assemble_vector!(b_free, name, lfuncs, dof, num_dependent_vars)
    end
    u_fix = assign_boundary_values(dof, num_dependent_vars)
    return b_free, u_fix
end

function validate_linear_functional(name::String, lfuncs::Vector{<:Tuple},
	dof::DegreesOfFreedom)
    mesh = dof.mesh
    physical_groups = mesh.gmodel.physical_groups
    if !(name in keys(physical_groups))
	error("$name in linear functional $m must be a physical name")
    end
    elt_type = mesh.elt_type_in[name]
    properties = mesh.elt_properties[elt_type]
    mesh_elt_name = properties.element_name
    for (m, lfunc) in enumerate(lfuncs)
	elt_vector!, source, params... = lfunc
        lfunc_elt_name, lfunc_elt_dof = elt_vector!()
        if lfunc_elt_name ≠ mesh_elt_name
	    error("""$(m)th linear functional for '$name' requires elements 
		  of type '$lfunc_elt_name' but '$name' has elements 
		  of type '$mesh_elt_name'""")
        end
	if lfunc_elt_dof ≠ dof.elt_dof[elt_type]
	    error("""$(m)th linear functional for '$name' assumes element 
	          degrees of freedom $lfunc_elt_dof but '$name' provides
	          $(dof.elt_dof[elt_type])""")
	end
    end
end

"""
    assemble_vector!(global_vec, name, lfuncs, dof, num_dependent_vars)

Assembles the contribution to the load vector from all linear functionals
from the vector `lfuncs` over the spatial domain `name`.
"""
function assemble_vector!(global_vec::Vector{Float64}, name::String, 
	lfuncs::Vector{<:Tuple}, dof::DegreesOfFreedom, 
	num_dependent_vars::Int64)

    each_elt_index, num_free, _, elt_global_dof, num_elt_dof, 
    coord = prepare_assembly(name, dof, num_dependent_vars)
    local_vec = Vector{Float64}(undef, num_elt_dof * num_dependent_vars)
    Σ_local_vec = similar(local_vec)
    for l in each_elt_index
	elt_node_coord!(coord, elt_global_dof, name, l, dof)
	fill!(Σ_local_vec, 0.0)
	for lfunc in lfuncs
	    elt_vector!, f, params... = lfunc
	    elt_vector!(local_vec, coord, f, params...)
	    Σ_local_vec .+= local_vec
	end
	for j in eachindex(elt_global_dof)
	    r = elt_global_dof[j]
	    if r > num_free 
		continue
	    end
            for p = 1:num_dependent_vars
		rp = r + (p-1)*num_free
		jp = j + (p-1)*num_elt_dof
		global_vec[rp] += Σ_local_vec[jp]
            end
	end
    end
end

function prepare_assembly(name::String, dof::DegreesOfFreedom,
                          num_dependent_vars::Int64)
    mesh = dof.mesh
    num_free, num_fixed = dof.num_free, dof.num_fixed
    elt_tags = mesh.elt_tags_in[name]
    each_elt_index = eachindex(elt_tags)
    elt_type = mesh.elt_type_in[name]
    elt_dof = dof.elt_dof[elt_type]
    num_elt_dof = lastindex(elt_dof)
    elt_prop = mesh.elt_properties[elt_type]
    num_elt_nodes = elt_prop.num_nodes
    coord = Matrix{Float64}(undef, (2, num_elt_nodes))
    elt_global_dof = Vector{Int64}(undef, num_elt_dof)
    return each_elt_index, num_free, num_fixed, elt_global_dof, 
           num_elt_dof, coord
end

function assign_boundary_values(dof::DegreesOfFreedom,
	num_dependent_vars::Int64)
    num_free, num_fixed = dof.num_free, dof.num_fixed
    node_tag_to_dof = dof.node_tag_to_dof
    mesh = dof.mesh
    u_fix = zeros(num_fixed * num_dependent_vars)
    valg = Vector{Float64}(undef, num_dependent_vars)
    for (name, g) in dof.essential_bc
        elt_node_tags = mesh.elt_node_tags_in[name]
        num_elts = size(elt_node_tags, 2)
        elt_type = mesh.elt_type_in[name]
        edof = dof.elt_dof[elt_type]
        if g isa Function
            for j = 1:num_elts
                for i in edof
                    nd = elt_node_tags[i,j]
                    k = node_tag_to_dof[nd] - num_free
                    x = mesh.coord[1,nd]
                    y = mesh.coord[2,nd]
                    valg .= g(x,y)
                    for p = 1:num_dependent_vars
                        u_fix[k+(p-1)*num_fixed] = valg[p]
                    end
                end
            end
        else # Constant boundary condition
            for j = 1:num_elts
                for i in edof
                    nd = elt_node_tags[i,j]
                    k = node_tag_to_dof[nd] - num_free
                    valg .= g
                    for p = 1:num_dependent_vars
                        u_fix[k+(p-1)*num_fixed] = valg[p]
                    end
                end
            end
        end
    end
    return u_fix
end

function assemble_matrix(dof::DegreesOfFreedom,
	bilinear_forms::Dict, num_dependent_vars=1)
    num_free, num_fixed = dof.num_free, dof.num_fixed
    A = sparse(Int64[], Int64[], Float64[], 
	       num_free*num_dependent_vars, 
	       (num_free+num_fixed)*num_dependent_vars)
    for (name, bforms) in bilinear_forms
	if bforms isa Tuple
	    bforms = [bforms]
	end
	validate_bilinear_form(name, bforms, dof)
	A += assemble_matrix(name, bforms, dof, num_dependent_vars)
    end
    A_free = A[:,1:num_free*num_dependent_vars]
    A_fix = A[:,num_free*num_dependent_vars+1:end]
    return A_free, A_fix
end

function assemble_matrix(name::String, bforms::Vector{<:Tuple}, 
	dof::DegreesOfFreedom, num_dependent_vars::Int64)
    each_elt_index, num_free, num_fixed, elt_global_dof, num_elt_dof, 
    coord = prepare_assembly(name, dof, num_dependent_vars)
    local_mat = Matrix{Float64}(undef, num_elt_dof*num_dependent_vars, 
				       num_elt_dof*num_dependent_vars)
    Σ_local_mat = similar(local_mat)
    num_elts = lastindex(each_elt_index)

    max_num_entries = num_dependent_vars^2 * num_elts * num_elt_dof^2
    I = Vector{Int64}(undef, max_num_entries)
    J = similar(I)
    V = Vector{Float64}(undef, max_num_entries)
    n = 0
    for l in each_elt_index
	elt_node_coord!(coord, elt_global_dof, name, l, dof)
	fill!(Σ_local_mat, 0.0)
	for bform in bforms
	    elt_matrix!, coef, params... = bform
            fill!(local_mat, 0.0)
	    elt_matrix!(local_mat, coord, coef, params...)
	    Σ_local_mat .+= local_mat
	end
	for i = eachindex(elt_global_dof)
            r = elt_global_dof[i]
	    if r > num_free
	        continue
	    end
            for p = 1:num_dependent_vars
	        rp = r + (p-1)*num_free
	        ip = i + (p-1)*num_elt_dof
		for j in eachindex(elt_global_dof)
                    s = elt_global_dof[j]
	            for q = 1:num_dependent_vars
		        if s ≤ num_free
		            offset = (q-1) * num_free
		        else
			    offset = ( (num_dependent_vars-1) * num_free 
				     + (q-1) * num_fixed )
		        end
		        sq = s + offset
		        jq = j + (q-1)*num_elt_dof
			n += 1
                        I[n] = rp
                        J[n] = sq
                        V[n] = Σ_local_mat[ip,jq]
		    end
	        end
	    end
        end
    end
    In = view(I, 1:n)
    Jn = view(J, 1:n)
    Vn = view(V, 1:n)
    return sparse(In, Jn, Vn, num_free*num_dependent_vars, 
		      (num_free+num_fixed)*num_dependent_vars)
end

function validate_bilinear_form(name::String, bforms::Vector{<:Tuple}, 
	dof::DegreesOfFreedom)
    mesh = dof.mesh
    physical_groups = mesh.gmodel.physical_groups
    if !(name in keys(physical_groups))
	error("$name in bilinear form $m must be a physical name")
    end
    elt_type = mesh.elt_type_in[name]
    properties = mesh.elt_properties[elt_type]
    mesh_elt_name = properties.element_name
    for (m, bform) in enumerate(bforms)
        elt_matrix!, source, params... = bform
        bform_elt_name, bform_elt_dof = elt_matrix!()
        if bform_elt_name ≠ mesh_elt_name
	    error("""$(m)th bilinear form for '$name' requires elements of 
		     type '$bform_elt_name' but '$name' has elements of 
		     type '$mesh_elt_name'""")
        end
	if bform_elt_dof ≠ dof.elt_dof[elt_type]
	    error("""$(m)th bilinear form for '$name' assumes element degrees
		     of freedom $bform_elt_dof but '$name' provides
		     $(dof.elt_dof[elt_type])""")
	end
    end
end

function average_field(uh::Vector{Float64}, name::String, dof::DegreesOfFreedom)
    each_elt_index, num_free, _, elt_global_dof, num_elt_dof,
    coord = prepare_assembly(name, dof, 1)
    ∫∫u = 0.0
    domain_area = 0.0
    for l in each_elt_index
        elt_node_coord!(coord, elt_global_dof, name, l, dof)
        Σ = 0.0
        for j in eachindex(elt_global_dof)
            r = elt_global_dof[j]
            Σ += uh[r]
        end
        A = triangle_area(coord)
        ∫∫u += A * Σ / 3
        domain_area += A
    end
    return ∫∫u, domain_area
end

function triangle_area(coord::Matrix{Float64})
    v1 = coord[1,1] - coord[1,3]
    v2 = coord[2,1] - coord[2,3]
    w1 = coord[1,2] - coord[1,3]
    w2 = coord[2,2] - coord[2,3]
    # overwrite B[1:2,1:2] with its inverse transpose
    area = abs(v1 * w2 - v2 * w1) / 2
end

end # module
