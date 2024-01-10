module MeshGen

import Gmsh: gmsh
import Base: getindex
import ..GeometryModel, ..DegreesOfFreedom, ..FEMesh, ..FuncOrConst,
       ..Edge, ..EdgeEnumeration, ..PRED
import ..Utils: diameter

function GeometryModel(name::String, verbosity=3)
    gmsh.initialize()
    gmsh.option.set_number("General.Verbosity", verbosity)
    stem, ext = splitext(name)
    if ext == ""
	filename = string(name, ".geo")
    else
	filename = name
    end
    if isfile(filename)
        gmsh.open(filename)
    else
        gmsh.finalize()
	error("File $filename not found")
    end
    entities = gmsh.model.get_entities()
    physical_groups = Dict{String,Tuple{Int32,Int32}}()
    entities_in = Dict{String,Vector{Int32}}()
    for (dim,tag) in gmsh.model.get_physical_groups()
	physical_name = gmsh.model.get_physical_name(dim, tag)
	physical_groups[physical_name] = (dim, tag)
	entities_in[physical_name] = (
            gmsh.model.get_entities_for_physical_group(dim,tag))
    end
    gmsh.finalize()
    return GeometryModel(filename, entities, physical_groups, entities_in)
end

function getindex(dof::DegreesOfFreedom, j::Int64)
    return dof.node_tag[j]
end

function DegreesOfFreedom(mesh::FEMesh,
	essential_bc::Vector{<:Tuple},
	elt_dof::Dict{Int32,Vector{Int64}}
	= default_elt_dof(mesh))
    # Determine the status of all nodes: status[nd] is
    #     0 if nd is not part of the physical mesh
    #     1 if nd is a not a free or fixed node (used for nonconforming elt)
    #     2 if nd is a free node
    #     3 if nd is a fixed node
    num_nodes = size(mesh.coord, 2)
    physical_groups = mesh.gmodel.physical_groups
    elt_type_in = mesh.elt_type_in
    elt_node_tags_in = mesh.elt_node_tags_in
    status = zeros(Int64, num_nodes)
    for (name, dim_tag) in physical_groups
	dim, tag = dim_tag
	if dim == 2
	    et = elt_type_in[name]
	    edof = elt_dof[et]
            elt_node_tags = elt_node_tags_in[name]
            num_elt_nodes, num_elts = size(elt_node_tags)
            ndof = Int64[]
            for i = 1:num_elt_nodes
		if !(i in edof)
                    push!(ndof, i)
                end
            end
	    for j = 1:num_elts
		for i in ndof
		    nd = elt_node_tags[i,j]
                    status[nd] = 1
		end
		for i in edof
		    nd = elt_node_tags[i,j]
                    status[nd] = 2
		end
	    end
	end
    end
    # Find every fixed node
    for (m,ebc) in enumerate(essential_bc)
	name, g = validate_essential_bc(ebc, m, mesh)
        elt_node_tags = elt_node_tags_in[name]
	num_elts = size(elt_node_tags, 2)
        et = elt_type_in[name]
	edof = elt_dof[et]
	for j = 1:num_elts
	    for i in edof
		nd = elt_node_tags[i,j]
		status[nd] = 3
	    end
	end
    end
    num_free  = 0
    num_fixed = 0
    for nd = 1:num_nodes
	if status[nd] == 2
	    num_free += 1
	elseif status[nd] == 3
	    num_fixed += 1
	end
    end
    num_phys_nodes = num_free + num_fixed
    # Enumerate the physical nodes with the free preceding the fixed
    node_tag = Vector{Int64}(undef, num_phys_nodes)
    k_free = 0
    k_fixed = num_free
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
    for k = 1:num_free
	nd = node_tag[k]
	node_tag_to_dof[nd] = k
    end
    for k = num_free+1:num_phys_nodes
	nd = node_tag[k]
	node_tag_to_dof[nd] = k
    end
    return DegreesOfFreedom(mesh, essential_bc, node_tag, num_free, 
			    num_fixed, node_tag_to_dof, elt_dof)
end

function FEMesh(gmodel::GeometryModel, hmax::Float64; 
	        order=1, refinements=0, verbosity=3, save_msh_file=false)
    gmsh.initialize()
    gmsh.option.set_number("General.Verbosity", verbosity)
    # Generate mesh and save to file
    filename = gmodel.name
    if isfile(filename)
        gmsh.open(filename)
    else
        gmsh.finalize()
	error("File $filename not found")
    end
    gmsh.option.set_number("Mesh.MeshSizeMax", hmax)
    gmsh.option.set_number("Mesh.ElementOrder", order)
    gmsh.model.mesh.generate(2)

    pathvec = splitpath(filename)
    name, ext = splitext(pathvec[end]) # strip everything except stem
    coord, elt_type_in, elt_tags_in, elt_node_tags_in, elt_properties =
        generate_mesh(gmodel, order, verbosity)
    if refinements==0
	mesh = FEMesh(gmodel, coord, elt_type_in, elt_tags_in,
		  elt_node_tags_in, elt_properties)
	if save_msh_file
	    write_msh_file(name, verbosity)
	end
    else
	mesh = Vector{FEMesh}(undef, refinements+1)
	mesh[1] = FEMesh(gmodel, coord, elt_type_in, elt_tags_in,
		  elt_node_tags_in, elt_properties)
	if save_msh_file
	    name_1 = string(name, 1)
	    write_msh_file(name_1, verbosity)
	end
	for l = 2:refinements+1
	    gmsh.model.mesh.refine()
	    gmsh.model.mesh.setOrder(order)
            coord, elt_type_in, elt_tags_in, elt_node_tags_in, elt_properties =
                generate_mesh(gmodel, order, verbosity)
	    mesh[l] = FEMesh(gmodel, coord, elt_type_in, elt_tags_in,
                             elt_node_tags_in, elt_properties)
	    if save_msh_file
	        name_l = string(name, l)
  	        write_msh_file(name_l, verbosity)
            end
	end
    end
    gmsh.finalize()
    return mesh
end

function write_msh_file(name::String, verbosity::Int64)
    filename = string(name, ".msh")
    if verbosity ≥ 3
        println("Writing mesh to ", filename)
    end
    gmsh.write(filename)
end

function generate_mesh(gmodel::GeometryModel, 
		       order::Int64, verbosity::Int64)
    physical_groups  = gmodel.physical_groups
    # Get node coordinates, dropping z=0
    node_tags, pre_coord, parametric_coord = gmsh.model.mesh.get_nodes()
    pre_coord = reshape(pre_coord, (3, lastindex(node_tags)))
    pre_coord = pre_coord[1:2,:]
    # Reorder to be consistent with elt_node_tags
    coord = similar(pre_coord)
    for j in eachindex(node_tags)
	nd = node_tags[j]
	coord[:,nd] = pre_coord[:,j]
    end
    
    elt_type_in      = Dict{String,Int32}()
    elt_tags_in      = Dict{String,Vector{Int64}}()
    elt_node_tags_in = Dict{String,Array{Int64}}()
    all_elt_types = Int32[]
    for (phys_name, dim_tag) in physical_groups
	dim, tag = dim_tag
	elt_type = Int32[]
	elt_tags = Int64[]
	elt_node_tags = Int64[]
	for etag in gmodel.entities_in[phys_name]
	    elt_data = gmsh.model.mesh.get_elements(dim, etag)
	    if length(elt_data[1]) > 1
		error("Each entity may have only a single element type")
	    end
	    push!(elt_type, elt_data[1][1])
	    append!(elt_tags, elt_data[2][1])
	    append!(elt_node_tags, elt_data[3][1])
	end
	et = elt_type[1]
	for j = 2:length(elt_type)
	    if et ≠ elt_type[j]
		error("Each physical group may have only a single element type")
	    end
	end
	elt_type_in[phys_name] = et
	elt_tags_in[phys_name] = elt_tags
	elt_node_tags_in[phys_name] = elt_node_tags
	push!(all_elt_types, et)
    end
    elt_properties = Dict{Int32,NamedTuple}()
    for et in Set(all_elt_types)
        props = gmsh.model.mesh.get_element_properties(et)
        elt_properties[et] = (element_name=props[1],
                              dim=props[2],
                              order=props[3],
                              num_nodes=props[4],
                              local_node_coord=props[5],
                              num_primary_nodes=props[6])
    end
    mat_elt_node_tags_in = Dict{String,Matrix{Int64}}()
    for (name, dim_tag) in physical_groups
	num_elts = lastindex(elt_tags_in[name])
        elt_type = elt_type_in[name]
        elt_prop = elt_properties[elt_type]
        nodes_per_elt = elt_prop.num_nodes
	mat_elt_node_tags_in[name] = reshape(elt_node_tags_in[name],
                                            (elt_prop.num_nodes, num_elts))
    end
    return coord, elt_type_in, elt_tags_in, 
           mat_elt_node_tags_in, elt_properties 
end

#function Base.display(gmodel::GeometryModel)
#    gmsh.initialize()
#    filename = gmodel.name
#    if isfile(filename)
#        gmsh.open(filename)
#    else
#	error("File $filename not found")
#    end
#    gmsh.option.set_number("Geometry.Points", 0)
#    gmsh.option.set_number("Geometry.LabelType", 4)
#    gmsh.option.set_number("Geometry.CurveLabels", 1)
#    gmsh.option.set_number("Geometry.Surfaces", 1)
#    gmsh.option.set_number("Geometry.SurfaceLabels", 1)
#    gmsh.fltk.run()
#    gmsh.finalize()
#end

function default_elt_dof(mesh::FEMesh)
    # By default, all nodes in an element are potential degrees of freedom.
    elt_potential_dof = Dict{Int32,Vector{Int64}}()
    for (k, ep) in mesh.elt_properties
        elt_potential_dof[k] = collect(1:ep.num_nodes)
    end
    return elt_potential_dof
end

function reorder_dof_steklov!(dof::DegreesOfFreedom, B_bilinear_forms::Dict)
    # Renumber the degrees of freedom so that the free nodes on the parts
    # of the boundary where the Steklov boundary condition applies
    # precede all of the other free nodes (and hence also the fixed nodes).
    mesh = dof.mesh
    num_free, num_fixed = dof.num_free, dof.num_fixed
    status = zeros(Int64, num_free)
    for (name, bforms) in B_bilinear_forms
	if bforms isa Tuple
            bforms = [bforms]
        end
	for bform in bforms
	    elt_matrix!, coef = bform
	    dim, tag = mesh.gmodel.physical_groups[name]
	    if dim ≠ 1
	        error("$name is not part of the boundary")
	    end
	    elt_node_tags = mesh.elt_node_tags_in[name]
	    elt_type = mesh.elt_type_in[name]
	    elt_dof = dof.elt_dof[elt_type]
	    num_elts = size(elt_node_tags, 2)
	    for l = 1:num_elts
	        for nd in elt_node_tags[elt_dof,l]
  	            j = dof.node_tag_to_dof[nd]
                    if j <= num_free
                        status[j] = 1
                    end
	        end
            end
	end
    end
    num_steklov = 0
    for j = 1:num_free
	if status[j] == 1
	    num_steklov += 1
	end
    end
    k_steklov = 0
    k_other = num_steklov
    new_node_tag = Vector{Int64}(undef, num_free)
    for j = 1:num_free
        if status[j] == 0
            k_other += 1
            new_node_tag[k_other] = dof[j]
        else
            k_steklov += 1
            new_node_tag[k_steklov] = dof[j]
        end
    end
    for j = 1:num_free
        nd = new_node_tag[j]
        dof.node_tag[j] = nd
        dof.node_tag_to_dof[nd] = j
    end
    return num_steklov
end

function elt_node_coord!(coord::Matrix{Float64}, 
	                 elt_global_dof::Vector{Int64}, name::String,
			 l::Int64, dof::DegreesOfFreedom)
    mesh = dof.mesh
    elt_node_tags = mesh.elt_node_tags_in[name]
    for j = 1:size(coord,2)
	nd = elt_node_tags[j,l]
	coord[1,j] = mesh.coord[1,nd]
	coord[2,j] = mesh.coord[2,nd]
    end
    elt_type = mesh.elt_type_in[name]
    elt_dof = dof.elt_dof[elt_type]
    for j in eachindex(elt_dof)
	nd = elt_node_tags[elt_dof[j],l]
	elt_global_dof[j] = dof.node_tag_to_dof[nd]
    end
end

function elt_node_coord!(coord::Matrix{Float64}, 
	                 elt_dof_tags::Vector{Int64}, mesh::FEMesh)
    for j in eachindex(elt_dof_tags)
	nd = elt_dof_tags[j]
	coord[:,j] = mesh.coord[:,nd]
    end
end

function get_nodal_values(u_func::Function, dof::DegreesOfFreedom)
    mesh = dof.mesh
    node_tag = dof.node_tag
    u_vec = Vector{Float64}(undef, lastindex(node_tag))
    for j in eachindex(node_tag)
	nd = node_tag[j]
	x = mesh.coord[1,nd]
	y = mesh.coord[2,nd]
	u_vec[j] = u_func(x, y)
    end
    return u_vec
end

function get_nodal_values(u_func::Function, dof::DegreesOfFreedom,
                          num_dependent_vars::Int64)
    mesh = dof.mesh
    node_tag = dof.node_tag
    u_mat = Matrix{Float64}(undef, (lastindex(node_tag), num_dependent_vars))
    for j in eachindex(node_tag)
	nd = node_tag[j]
	x = mesh.coord[1,nd]
	y = mesh.coord[2,nd]
	u_mat[j,:] .= u_func(x, y)
    end
    return u_mat
end

function validate_essential_bc(ebc::Tuple, m::Int64, mesh::FEMesh)
    if length(ebc) != 2
	error("Essential BC $m must be a Tuple of length 2")
    end
    name, g = ebc
    physical_groups = mesh.gmodel.physical_groups
    if !(name in keys(physical_groups))
	error("$name in Essential BC $m must be a physical name")
    end
    dim, tag = physical_groups[name]
    if dim != 1
	error("$name in Essential BC $m must be part of the boundary")
    end
    if !(typeof(g) <: FuncOrConst || typeof(g) <: FuncOrConst)
	error("Second member of Essential BC $m must be a function "*
	      "or Float64 constant")
    end
    return name, g
end

function max_elt_diameter(mesh::FEMesh)
    coord = Matrix{Float64}(undef, (2,3))
    hmax = 0.0
    for (name, dim_tag) in mesh.gmodel.physical_groups
        dim, tag = dim_tag
        if dim == 2
            elt_tags = mesh.elt_tags_in[name]
            elt_node_tags = mesh.elt_node_tags_in[name]
            for l in eachindex(elt_tags)
                elt_node_coord!(coord, elt_node_tags[1:3,l], mesh)
                d = diameter(coord)
                hmax = max(d, hmax)
            end
        end
    end
    return hmax
end

function EdgeEnumeration(mesh::FEMesh)
    num_nodes = size(mesh.coord, 2)
    isnew = fill(true, num_nodes)
    midpt_to_edge = fill(-1, num_nodes)
    edge = Edge[]
    i = 0
    for (name, dim_tag) in mesh.gmodel.physical_groups
	dim, tag = dim_tag
	if dim == 2
	    elt_tags = mesh.elt_tags_in[name]
	    elt_node_tags = mesh.elt_node_tags_in[name]
	    elt_type = mesh.elt_type_in[name]
	    elt_properties = mesh.elt_properties[elt_type]
	    @assert elt_properties.element_name == "Triangle 6"
	    for l in eachindex(elt_tags)
		for k = 1:3
		    midpt = elt_node_tags[k+3,l]
		    if isnew[midpt]
			isnew[midpt] = false
			i += 1
		        midpt_to_edge[midpt] = i
			left_phys_grp = tag
			left_triangle = l
			left_opposite = PRED[k]
			next_edge = Edge(midpt, left_phys_grp, left_triangle,
					 left_opposite, -1, -1, -1)
			push!(edge, next_edge)
		    else
			j = midpt_to_edge[midpt]
			edge[j].right_phys_grp = tag
			edge[j].right_triangle = l
			edge[j].right_opposite = PRED[k]
		    end
		end
	    end
	end
    end
    return EdgeEnumeration(edge, midpt_to_edge)
end

end # module
