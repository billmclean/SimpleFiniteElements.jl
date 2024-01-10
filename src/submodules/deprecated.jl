function assemble_vector(dof::DegreesOfFreedom,
        linear_functionals::Vector{<:Tuple}, num_dependent_vars::Int64=1)
    mesh = dof.mesh
    num_free, node_tag_to_dof = dof.num_free, dof.node_tag_to_dof
    # Accumulate the global vector in b_free
    b_free = zeros(num_free * num_dependent_vars)
    for (m, lfunc)  in enumerate(linear_functionals)
        name, elt_vector!, source, params... =
            validate_linear_functional(lfunc, m, mesh)
        # Contribution from current functional
        assemble_vector!(b_free, name, elt_vector!, source, dof,
                         num_dependent_vars, params...)
    end
    # Now assign boundary values at the fixed nodes
    num_fixed = dof.num_fixed
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
    return b_free, u_fix
end

function validate_bilinear_form(bform::Tuple, m::Int64, mesh::FEMesh)
    if length(bform) < 3
        error("Bilinear form $m must be a tuple of length 3 or more")
    end
    name, elt_matrix!, source, params... = bform
    physical_groups = mesh.gmodel.physical_groups
    if !(name in keys(physical_groups))
        error("$name in bilinear form $m must be a physical name")
    end
    elt_type = mesh.elt_type_in[name]
    properties = mesh.elt_properties[elt_type]
    bform_elt_name = elt_matrix!()
    mesh_elt_name = properties.element_name
    if bform_elt_name ≠ mesh_elt_name
        error("Bilinear form $m requires elements of type '$bform_elt_name'"
              * " but '$name' has elements of type '$mesh_elt_name'")
    end
    return bform
end

function validate_linear_functional(lfunc::Tuple, m::Int64, mesh::FEMesh)

    if length(lfunc) < 3
        error("Linear functional $m must be a tuple of length 3 or more")
    end
    name, elt_vector!, source, params... = lfunc
    elt_type = mesh.elt_type_in[name]
    properties = mesh.elt_properties[elt_type]
    lfunc_elt_name = elt_vector!()
    mesh_elt_name = properties.element_name
    if lfunc_elt_name ≠ mesh_elt_name
        error("Linear functional $m requires elements of type '$lfunc_elt_name'"
              * " but '$name' has elements of type '$mesh_elt_name'")
    end
    return lfunc
end

function assemble_vector!(global_vec::Vector{Float64},
                          name::String, elt_vector!::Function,
                          f::FuncOrConst, dof::DegreesOfFreedom,
                          num_dependent_vars::Int64, params...)

    each_elt_index, num_free, _, elt_global_dof, num_elt_dof,
    coord = prepare_assembly(name, dof, num_dependent_vars)
    local_vec = Vector{Float64}(undef, num_elt_dof * num_dependent_vars)
    for l in each_elt_index
        elt_node_coord!(coord, elt_global_dof, name, l, dof)
        elt_vector!(local_vec, coord, f, params...)
        for j in eachindex(elt_global_dof)
            r = elt_global_dof[j]
            if r > num_free
                continue
            end
            for p = 1:num_dependent_vars
                rp = r + (p-1)*num_free
                jp = j + (p-1)*num_elt_dof
                global_vec[rp] += local_vec[jp]
            end
        end
    end
end

function assemble_matrix(dof::DegreesOfFreedom,
        bilinear_forms::Vector, num_dependent_vars=1)
    mesh = dof.mesh
    num_free, num_fixed = dof.num_free, dof.num_fixed
    A = sparse(Int64[], Int64[], Float64[],
               num_free*num_dependent_vars,
               (num_free+num_fixed)*num_dependent_vars)
    for (m, bform) in enumerate(bilinear_forms)
        name, elt_matrix!, coef, params... =
            validate_bilinear_form(bform, m, mesh)
        A += assemble_matrix(name, elt_matrix!, coef, dof,
                             num_dependent_vars, params...)
    end
    A_free = A[:,1:num_free*num_dependent_vars]
    A_fix = A[:,num_free*num_dependent_vars+1:end]
    return A_free, A_fix
end

function assemble_matrix(name::String, elt_matrix!::Function,
                         coef::FuncOrConst, dof::DegreesOfFreedom,
                         num_dependent_vars, params...)
    each_elt_index, num_free, num_fixed, elt_global_dof, num_elt_dof,
    coord = prepare_assembly(name, dof, num_dependent_vars)
    local_mat = Matrix{Float64}(undef, num_elt_dof*num_dependent_vars,
                                       num_elt_dof*num_dependent_vars)
    num_elts = lastindex(each_elt_index)

    max_num_entries = num_dependent_vars^2 * num_elts * num_elt_dof^2
    I = Vector{Int64}(undef, max_num_entries)
    J = similar(I)
    V = Vector{Float64}(undef, max_num_entries)
    n = 0
    for l in each_elt_index
        elt_node_coord!(coord, elt_global_dof, name, l, dof)
        elt_matrix!(local_mat, coord, coef, params...)
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
                        V[n] = local_mat[ip,jq]
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


