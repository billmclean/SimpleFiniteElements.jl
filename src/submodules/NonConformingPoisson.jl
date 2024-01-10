module NonConformingPoisson

import ..GeometryModel, ..FEMesh, ..DegreesOfFreedom, ..SUCC, ..PRED
import ..Utils: barycentric, centroid, quadrature_points, nc_shape_func!
import ..Utils: barycentric!
import ..FEM: prepare_assembly
import ..MeshGen: elt_node_coord!
import LinearAlgebra: lmul!

const ELT_DOF =  Dict(
    Int32(9) => [ 4, 5, 6 ],
    Int32(8) => [ 3 ])

function ∫∫a_∇u_dot_∇v!(A::Matrix{Float64}, coord::Matrix{Float64},
                      coef::Float64)
    b, area = barycentric(coord[:,1:3])
    C = 4 * area * coef
    for k = 1:3
	pk = PRED[k]
	for j = k:3
	    pj = PRED[j]
	    A[j,k] = C * ( b[1,pk] * b[1,pj] + b[2,pk] * b[2,pj] )
	end
    end
    for k = 2:3
	for j = 1:k-1
	    A[j,k] = A[k,j]
	end
    end
end

function ∫∫a_∇u_dot_∇v!(A::Matrix{Float64}, coord::Matrix{Float64},
                      coef::Function, params...)
    b, area = barycentric(coord[:,1:3])
    Σ = 0.0
    for p = 4:6
	x = coord[1,p]
	y = coord[2,p]
	Σ += coef(x, y, params...)
    end
    Σ *= (4/3) * area
    for k = 1:3
	pk = PRED[k]
	for j = k:3
	    pj = PRED[j]
	    A[j,k] = Σ * ( b[1,pk] * b[1,pj] + b[2,pk] * b[2,pj] )
	end
    end
    for k = 2:3
	for j = 1:k-1
	    A[j,k] = A[k,j]
	end
    end
end

∫∫a_∇u_dot_∇v!() = ("Triangle 6", [4, 5, 6])

function ∫∫c_u_v!(M::Matrix{Float64}, coord::Matrix{Float64}, coef::Float64)
    b, area = barycentric(coord[:,1:3])
    fill!(M, 0.0)
    for p = 1:3
	M[p,p] = coef * area/3
    end
end

function ∫∫c_u_v!(M::Matrix{Float64}, coord::Matrix{Float64}, 
	coef::Function, params...)
    b, area = barycentric(coord[:,1:3])
    fill!(M, 0.0)
    for p = 1:3
	k = p + 3
	x, y = coef[1,k], coef[2,k]
	M[p,p] = (coef * area/3) * f(x, y, params...)
    end
end

∫∫c_u_v!() = ("Triangle 6", [4, 5, 6])

function ∫∫f_v!(v::Vector{Float64}, coord::Matrix{Float64}, f::Float64)
    b, area = barycentric(coord)
    fill!(v, f * area / 3)
end

function ∫∫f_v!(v::Vector{Float64}, coord::Matrix{Float64}, 
	f::Function, params...)
    b, area = barycentric(coord)
    for k = 1:3
	p = k + 3
	x = coord[1,p]
	y = coord[2,p]
	v[k] = (area/3) * f(x, y, params...)
    end
end

∫∫f_v!() = ("Triangle 6", [4, 5, 6])

function ∫c_u_v!(M::Matrix{Float64}, coord::Matrix{Float64}, coef::Float64)
    x1, y1 = coord[1,1], coord[2,1]
    x2, y2 = coord[1,2], coord[2,2]
    lenE = hypot(x2-x1, y2-y1)
    M[1,1] = lenE * coef
end

function ∫c_u_v!(M::Matrix{Float64}, coord::Matrix{Float64}, 
	coef::Function, params...)
    x1, y1 = coord[1,1], coord[2,1]
    x2, y2 = coord[1,2], coord[2,2]
    x3, y3 = coord[1,2], coord[2,2]
    lenE = hypot(x2-x1, y2-y1)
    M[1,1] = lenE * coef(x3, y3, params...)
end

∫c_u_v!() = ("Line 3", [3])

function ∫g_v!(v::Vector{Float64}, coord::Matrix{Float64}, g::Float64)
    x1, y1 = coord[1,1], coord[2,1]
    x2, y2 = coord[1,2], coord[2,2]
    lenE = hypot(x2-x1, y2-y1)
    v[1] = lenE * g
end

function ∫g_v!(v::Vector{Float64}, coord::Matrix{Float64}, 
	g::Function, params...)
    x1, y1 = coord[1,1], coord[2,1]
    x2, y2 = coord[1,2], coord[2,2]
    x3, y3 = coord[1,3], coord[2,3]
    lenE = hypot(x2-x1, y2-y1)
    v[1] = lenE * g(x3, y3, params...)
end

∫g_v!() = ("Line 3", [3])

function gmsh2pyplot(dof::DegreesOfFreedom, 
	uh_midpt::Vector{Float64})
    mesh = dof.mesh
    physical_groups = mesh.gmodel.physical_groups
    coord = mesh.coord
#    num_mesh_nodes = length(dof.node_tag)
    num_triangles = 0
    for (name, dim_tag) in physical_groups
        dim, tag = dim_tag
        if dim == 2
            et = mesh.elt_type_in[name]
            @assert mesh.elt_properties[et].element_name == "Triangle 6"
	    elt_tags = mesh.elt_tags_in[name]
	    num_triangles += length(elt_tags)
        end
    end
    num_new_nodes = 3 * num_triangles
    x = Vector{Float64}(undef, num_new_nodes)
    y = similar(x)
    uh = similar(x)
    triangles = Matrix{Float64}(undef, (num_triangles,3))
    i = k = 0
    for (name, dim_tag) in physical_groups
        dim, tag = dim_tag
        if dim == 2
	    elt_node_tags = mesh.elt_node_tags_in[name]
	    for p = 1:size(elt_node_tags, 2)
		vertex_nd = elt_node_tags[1:3,p]
		midpt_nd = elt_node_tags[4:6,p]
		midpt_dof = dof.node_tag_to_dof[midpt_nd]
		i += 1
		for j = 1:3
		    k += 1
		    nd = vertex_nd[j]
		    x[k] = coord[1,nd]
		    y[k] = coord[2,nd]
		    triangles[i,j] = k - 1
		    uh[k] = ( uh_midpt[midpt_dof[j]] 
			    - uh_midpt[midpt_dof[SUCC[j]]] 
			    + uh_midpt[midpt_dof[PRED[j]]] )
		end
	    end
	end
    end
    return x, y, triangles, uh
end

function error_norms(uh::Vector{Float64}, u::Function, ∇u::Function,
                     dof::DegreesOfFreedom, level::Int64, params...)
    mesh = dof.mesh
    L2err, H1err = 0.0, 0.0
    ψ = Vector{Float64}(undef, 3)
    b = Matrix{Float64}(undef, 2, 3)
    cntrd = Vector{Float64}(undef, 2)
    ∇uh = similar(cntrd)
    for (name, elt_type) in mesh.elt_type_in
        # Ignore unless element type is "Triangle 6".
        if elt_type ≠ 9
	    continue
	end
        each_elt_index, _, _, elt_global_dof, _,
        coord = prepare_assembly(name, dof, 1)
        for l in each_elt_index
            elt_node_coord!(coord, elt_global_dof, name, l, dof)
            area = barycentric!(b, cntrd, coord)
            x, y = quadrature_points(level, coord)
            fill!(∇uh, 0.0) # ∇uh is constant on each element
            for j = eachindex(elt_global_dof)
                r = elt_global_dof[j]
		j₋ = PRED[j]
		∇uh[1] += uh[r] * (-2 * b[1,j₋])
		∇uh[2] += uh[r] * (-2 * b[2,j₋])
            end
            Σ_L2, Σ_H1 = 0.0, 0.0
            for k in eachindex(x)
                nc_shape_func!(ψ, x[k], y[k], b, cntrd)
                uh_k = 0.0
                for j = eachindex(elt_global_dof)
                    r = elt_global_dof[j]
                    uh_k += uh[r] * ψ[j]
                end
                Σ_L2 +=  (uh_k - u(x[k], y[k], params...))^2
                ∇u_k = ∇u(x[k], y[k], params...)
                Σ_H1 += (∇uh[1] - ∇u_k[1])^2 + (∇uh[2] - ∇u_k[2])^2
            end
            L2err += (area / lastindex(x)) * Σ_L2
            H1err += (area / lastindex(x)) * Σ_H1
        end
    end
    L2err = sqrt(L2err)
    H1err = sqrt(H1err)
    return L2err, H1err
end

end # module
