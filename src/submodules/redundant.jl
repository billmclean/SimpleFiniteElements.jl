function FEMesh_midpoint_dof(gmodel::GeometryModel, hmax::Float64,
                essential_bc::Vector{EssentialBC}, verbosity=3)
    order = 2
    coord, elt_type_in, elt_tags_in, elt_node_tags_in, elt_properties =
        generate_mesh(gmodel, hmax, essential_bc, order, verbosity)
    num_nodes = size(coord, 2)
    physical_groups = gmodel.physical_groups
    node_tag, num_free_midpt, num_fixed_midpt, node_tag_to_dof, which_bc =
        assign_midpt_degrees_of_freedom(essential_bc, num_nodes, 
					physical_groups, elt_type_in,
					elt_node_tags_in, elt_properties)
    which_nodes = [4, 5, 6]
    dof = DegreesOfFreedom(node_tag, num_free_midpt, num_fixed_midpt,
                           node_tag_to_dof, which_bc, which_nodes)
    return FEMesh(gmodel, coord, elt_type_in, elt_tags_in,
                  elt_node_tags_in, elt_properties, dof, essential_bc)
end

function assign_midpt_degrees_of_freedom(essential_bc::Vector{EssentialBC},
        num_nodes::Int64, physical_groups::Dict{String,Tuple{Int32,Int32}},
	elt_type_in::Dict{String,Int32}, 
        elt_node_tags_in::Dict{String,Matrix{Int64}},
	elt_properties::Dict{Int32,NamedTuple})
    # Determine status of all nodes: status[nd] is
    #    0 if nd is not part of the physical mesh
    #    1 if nd is a vertex of the physical mesh
    #    2 if nd is a free midpoint of the physical mesh
    #    3 if nd is a fixed midpoint of the physical mesh
    status = zeros(Int64, num_nodes)
    for (name, dim_tag) in physical_groups
	dim, tag = dim_tag
	if dim == 2
	    elt_type = elt_type_in[name]
	    elt_prop = elt_properties[elt_type]
	    if elt_prop.element_name ≠ "Triangle 6"
		error("Elements in $name must be 6 node triangles")
	    end
	    elt_node_tags = elt_node_tags_in[name]
	    num_elts = size(elt_node_tags, 2)
#	    num_elts = div(length(elt_node_tags_in[name]), 6)
#	    ent = reshape(elt_node_tags_in[name], (6, num_elts))
	    for j = 1:num_elts
		for i = 1:3
		    nd = elt_node_tags[i,j]
		    status[nd] = 1 
		end
		for i = 4:6
		    nd = elt_node_tags[i,j]
		    status[nd] = 2 
		end
	    end
	end
    end
    for ebc in essential_bc
	name = ebc.physical_name
	elt_type = elt_type_in[name]
	elt_prop = elt_properties[elt_type]
	if elt_prop.element_name ≠ "Line 3"
            error("Elements in $name must be 3 node lines")
        end
        elt_node_tags = elt_node_tags_in[name]
        num_elts = size(elt_node_tags, 2)
	for j = 1:num_elts
	    nd = elt_node_tags[3,j]
	    status[nd] = 3
	end
    end
    # Enumerate the midpoint nodes with the free preceding the fixed
    num_free_midpt = 0
    num_fixed_midpt = 0
    num_vertex = 0
    for nd = 1:num_nodes
	if status[nd] == 1
	    num_vertex += 1
	elseif status[nd] == 2
	    num_free_midpt += 1
	elseif status[nd] == 3
	    num_fixed_midpt += 1
	end
    end
    num_midpt = num_free_midpt + num_fixed_midpt
    node_tag = Vector{Int64}(undef, num_midpt)
    k_free = 0
    k_fixed = num_free_midpt
    for nd = 1:num_nodes
	if status[nd] == 2
	    k_free += 1
	    node_tag[k_free] = nd
	elseif status[nd] == 3
	    k_fixed += 1
	    node_tag[k_fixed] = nd
	end
    end
    node_tag_to_dof = zeros(Int64, num_nodes)
    for k = 1:num_free_midpt
	nd = node_tag[k]
	node_tag_to_dof[nd] = k
    end
    for k = num_free_midpt+1:num_midpt
	nd = node_tag[k]
	node_tag_to_dof[nd] = k
    end
    which_bc = Vector{Int64}(undef, num_fixed_midpt)
    for j = 1:length(essential_bc)
	name = essential_bc[j].physical_name
	elt_node_tags = elt_node_tags_in[name]
        num_elts = size(elt_node_tags, 2)
	for l = 1:num_elts
	    nd = elt_node_tags[3,l]
	    k = node_tag_to_dof[nd] - num_free_midpt
	    which_bc[k] = j
	end
    end
    return node_tag, num_free_midpt, num_fixed_midpt, node_tag_to_dof, which_bc
end


function linear_system(dbvp::DiscreteBVP)
    mesh, dof = dbvp.mesh, dbvp.dof
    A = sparse(Int64[], Int64[], Float64[], 
	       dof.num_free, dof.num_free+dof.num_fixed)
    for bform in dbvp.bilinear_form
	A += assemble_matrix(bform.physical_name, bform.elt_matrix!, 
			     bform.coef, mesh, dof)
    end
    A_free = A[:,1:dof.num_free]
    A_fix = A[:,dof.num_free+1:end]
    b = zeros(dof.num_free)
    for lfunc  in dbvp.linear_functional
        assemble_vector!(b, lfunc.physical_name, lfunc.elt_vector!, 
			 lfunc.source, mesh, dof)
    end
    b -= A_fix * dbvp.u_fixed
    return A_free, b
end

function nodal_solution(dbvp::DiscreteBVP, u_free::Vector{Float64})
    dof = dbvp.dof
    num_free, num_fixed = dof.num_free, dof.num_fixed
    num_nodes = num_free + num_fixed
    uh = Vector{Float64}(undef, num_nodes)
    for j = 1:num_free
	nd = dof.node_tag[j]
	uh[nd] = u_free[j]
    end
    for j = 1:num_fixed
	nd = dof.node_tag[num_free+j]
	uh[nd] = dbvp.u_fixed[j]
    end
    return uh
end

function alt_nodal_solution(dbvp::DiscreteBVP, u_free::Vector{Float64})
    uh = append!(u_free, dbvp.u_fixed)
end

function matrix_pencil(dep::DiscreteEigenproblem)
    mesh = dep.mesh
    dof = dep.dof
    A = sparse(Int64[], Int64[], Float64[], 
	       dof.num_free, dof.num_free+dof.num_fixed)
    for bform in dep.A_bilinear_form
	A += assemble_matrix(bform.physical_name, bform.elt_matrix!, 
			     bform.coef, mesh, dof)
    end
    A_free = A[:,1:dof.num_free]
    B = sparse(Int64[], Int64[], Float64[], 
	       dof.num_free, dof.num_free+dof.num_fixed)
    for bform in dep.B_bilinear_form
	B += assemble_matrix(bform.physical_name, bform.elt_matrix!, 
			     bform.coef, mesh, dof)
    end
    B_free = B[:,1:dof.num_free]
    return A_free, B_free
end

function matrix_pencil(dsep::DiscreteSteklovEigenproblem; block_matrices=false)
    mesh = dsep.mesh
    dof = mesh.dof
    num_free = dof.num_free
    num_steklov = dsep.num_steklov
    A = sparse(Int64[], Int64[], Float64[], 
	       dof.num_free, dof.num_free+dof.num_fixed)
    for bform in dsep.A_bilinear_form
	A += assemble_matrix(bform.physical_name, bform.elt_matrix!, 
			     bform.coef, mesh)
    end
    B = sparse(Int64[], Int64[], Float64[], 
	       dof.num_free, dof.num_free+dof.num_fixed)
    for bform in dsep.B_bilinear_form
	B += assemble_matrix(bform.physical_name, bform.elt_matrix!, 
			     bform.coef, mesh)
    end
    if block_matrices
        A_free = Matrix{SparseMatrixCSC}(undef, 2, 2)
        A_free[1,1] = A[1:num_steklov,1:num_steklov]
        A_free[2,1] = A[num_steklov+1:num_free,1:num_steklov]
        A_free[1,2] = A[1:num_steklov,num_steklov+1:num_free]
        A_free[2,2] = A[num_steklov+1:num_free,num_steklov+1:num_free]
        B_free = Matrix{SparseMatrixCSC}(undef, 2, 2)
        B_free[1,1] = B[1:num_steklov,1:num_steklov]
        B_free[2,1] = B[num_steklov+1:num_free,1:num_steklov]
        B_free[1,2] = B[1:num_steklov,num_steklov+1:num_free]
        B_free[2,2] = B[num_steklov+1:num_free,num_steklov+1:num_free]
        if nnz(B_free[1,2]) > 0 || nnz(B_free[2,1]) > 0 || nnz(B_free[2,2]) > 0
            error("Block matrix B should have non-zeros only in top left block")
        end
    else
        A_free = A[:,1:num_free]
	B_free = B[:,1:num_free]
    end
    return A_free, B_free
end

function nodal_solution(dep::AbstractEigenproblem, ϕ_free::Vector{Float64})
    dof = dep.dof
    num_free, num_fixed = dof.num_free, dof.num_fixed
    num_nodes = num_free + num_fixed
    ϕh = zeros(num_nodes)
    for j = 1:num_free
	nd = dof.node_tag[j]
	ϕh[nd] = ϕ_free[j]
    end
    return ϕh
end


function linear_system(dbvp::DiscreteBVP)
    mesh, dof = dbvp.mesh, dbvp.dof
    A = sparse(Int64[], Int64[], Float64[], 
	       dof.num_free, dof.num_free+dof.num_fixed)
    for bform in dbvp.bilinear_form
	A += assemble_matrix(bform.physical_name, bform.elt_matrix!, 
			     bform.coef, mesh, dof)
    end
    A_free = A[:,1:dof.num_free]
    A_fix = A[:,dof.num_free+1:end]
    b = zeros(dof.num_free)
    for lfunc  in dbvp.linear_functional
        assemble_vector!(b, lfunc.physical_name, lfunc.elt_vector!, 
			 lfunc.source, mesh, dof)
    end
    b -= A_fix * dbvp.u_fixed
    return A_free, b
end

function nodal_solution(dbvp::DiscreteBVP, u_free::Vector{Float64})
    dof = dbvp.dof
    num_free, num_fixed = dof.num_free, dof.num_fixed
    num_nodes = num_free + num_fixed
    uh = Vector{Float64}(undef, num_nodes)
    for j = 1:num_free
	nd = dof.node_tag[j]
	uh[nd] = u_free[j]
    end
    for j = 1:num_fixed
	nd = dof.node_tag[num_free+j]
	uh[nd] = dbvp.u_fixed[j]
    end
    return uh
end

function alt_nodal_solution(dbvp::DiscreteBVP, u_free::Vector{Float64})
    uh = append!(u_free, dbvp.u_fixed)
end

function matrix_pencil(dep::DiscreteEigenproblem)
    mesh = dep.mesh
    dof = dep.dof
    A = sparse(Int64[], Int64[], Float64[], 
	       dof.num_free, dof.num_free+dof.num_fixed)
    for bform in dep.A_bilinear_form
	A += assemble_matrix(bform.physical_name, bform.elt_matrix!, 
			     bform.coef, mesh, dof)
    end
    A_free = A[:,1:dof.num_free]
    B = sparse(Int64[], Int64[], Float64[], 
	       dof.num_free, dof.num_free+dof.num_fixed)
    for bform in dep.B_bilinear_form
	B += assemble_matrix(bform.physical_name, bform.elt_matrix!, 
			     bform.coef, mesh, dof)
    end
    B_free = B[:,1:dof.num_free]
    return A_free, B_free
end

function matrix_pencil(dsep::DiscreteSteklovEigenproblem; block_matrices=false)
    mesh = dsep.mesh
    dof = mesh.dof
    num_free = dof.num_free
    num_steklov = dsep.num_steklov
    A = sparse(Int64[], Int64[], Float64[], 
	       dof.num_free, dof.num_free+dof.num_fixed)
    for bform in dsep.A_bilinear_form
	A += assemble_matrix(bform.physical_name, bform.elt_matrix!, 
			     bform.coef, mesh)
    end
    B = sparse(Int64[], Int64[], Float64[], 
	       dof.num_free, dof.num_free+dof.num_fixed)
    for bform in dsep.B_bilinear_form
	B += assemble_matrix(bform.physical_name, bform.elt_matrix!, 
			     bform.coef, mesh)
    end
    if block_matrices
        A_free = Matrix{SparseMatrixCSC}(undef, 2, 2)
        A_free[1,1] = A[1:num_steklov,1:num_steklov]
        A_free[2,1] = A[num_steklov+1:num_free,1:num_steklov]
        A_free[1,2] = A[1:num_steklov,num_steklov+1:num_free]
        A_free[2,2] = A[num_steklov+1:num_free,num_steklov+1:num_free]
        B_free = Matrix{SparseMatrixCSC}(undef, 2, 2)
        B_free[1,1] = B[1:num_steklov,1:num_steklov]
        B_free[2,1] = B[num_steklov+1:num_free,1:num_steklov]
        B_free[1,2] = B[1:num_steklov,num_steklov+1:num_free]
        B_free[2,2] = B[num_steklov+1:num_free,num_steklov+1:num_free]
        if nnz(B_free[1,2]) > 0 || nnz(B_free[2,1]) > 0 || nnz(B_free[2,2]) > 0
            error("Block matrix B should have non-zeros only in top left block")
        end
    else
        A_free = A[:,1:num_free]
	B_free = B[:,1:num_free]
    end
    return A_free, B_free
end

function nodal_solution(dep::AbstractEigenproblem, ϕ_free::Vector{Float64})
    dof = dep.dof
    num_free, num_fixed = dof.num_free, dof.num_fixed
    num_nodes = num_free + num_fixed
    ϕh = zeros(num_nodes)
    for j = 1:num_free
	nd = dof.node_tag[j]
	ϕh[nd] = ϕ_free[j]
    end
    return ϕh
end
