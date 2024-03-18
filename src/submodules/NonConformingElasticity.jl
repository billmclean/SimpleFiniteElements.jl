module NonConformingElasticity

import ..PRED, ..SUCC, ..FEMesh, ..DegreesOfFreedom, ..EdgeEnumeration
import ..Elasticity: elasticity_soln
import ..Utils: barycentric, centroid, quadrature_points, nc_shape_func!
import ..Utils: barycentric!
import ..MeshGen: elt_node_coord!
import ..NonConformingPoisson: ELT_DOF
import ..FEM: prepare_assembly
using PyPlot
import StaticArrays: SVector, SA

function visualise_soln(dof::DegreesOfFreedom, edenum::EdgeEnumeration,
        u1h::Vector{Float64}, u2h::Vector{Float64}, scale::Float64, figno=1)
    mesh = dof.mesh
    num_midpts = length(dof.node_tag)
    x = Vector{Float64}(undef, num_midpts)
    y = similar(x)
    for j = 1:num_midpts
	nd = dof.node_tag[j]
	x[j] = mesh.coord[1,nd]
	y[j] = mesh.coord[2,nd]
    end
    fig1 = figure(figno)
    if isnan(scale)
        quiver(x, y, u1h, u2h, angles="xy", scale_units="xy") # default scale
    else
        quiver(x, y, u1h, u2h, angles="xy", scale_units="xy", scale=scale)
    end
    xlabel("x")
    ylabel("y")
    axis("equal")

    gmodel = mesh.gmodel
    coord = mesh.coord
    node_tag_to_dof = dof.node_tag_to_dof

    triangle = Dict{Int32,Matrix{Int64}}()
    for (name, dim_tag) in gmodel.physical_groups
	dim, tag = dim_tag
	if dim == 2
	    triangle[tag] = mesh.elt_node_tags_in[name]
	end
    end

    function add_points!(x, y, x_ref, y_ref, edge)
	elt_node_tags = triangle[edge.left_phys_grp]
	p = edge.left_triangle
	k = edge.left_opposite
	midpt_nd = elt_node_tags[4:6,p]
	midpt_dof = dof.node_tag_to_dof[midpt_nd]
	u1h_vtx = vertex_values(u1h, midpt_dof)
	u2h_vtx = vertex_values(u2h, midpt_dof)
	for j in (SUCC[k], PRED[k])
  	    nd = elt_node_tags[j,p] 
            push!(x_ref, coord[1,nd])
	    push!(y_ref, coord[2,nd])
	    push!(x, coord[1,nd] + u1h_vtx[j]/scale)
	    push!(y, coord[2,nd] + u2h_vtx[j]/scale)
	end
	push!(x, NaN)
	push!(y, NaN)
    end

    figure(figno+1)
    for (name, dim_tag) in gmodel.physical_groups
        dim, tag = dim_tag
        x_ref = Float64[]
        y_ref = Float64[]
        x = Float64[]
        y = Float64[]
        if dim == 1
            elt_node_tags = mesh.elt_node_tags_in[name]
            num_elts = size(elt_node_tags, 2)
            for j = 1:num_elts
		midpt = elt_node_tags[3,j]
		i = edenum.midpt_to_edge[midpt]
		edge = edenum.edge[i]
		add_points!(x, y, x_ref, y_ref, edge)
            end
            plot(x_ref, y_ref, "k-", linewidth=0.5)
            plot(x, y, "b-")
            figure(figno)
            plot(x_ref, y_ref, "k-", linewidth=0.5)
            figure(figno+1)
        end
    end
    xlabel("x")
    ylabel("y")
    axis("equal")
end

function vertex_values(uh::Vector{Float64}, midpt_dof::Vector{Int64})
    uh_vtx = Vector{Float64}(undef, 3)
    for j = 1:3
        uh_vtx[j] = ( uh[midpt_dof[j]]
                    + uh[midpt_dof[PRED[j]]]
                    - uh[midpt_dof[SUCC[j]]] )
    end
    return uh_vtx
end

function ∫∫a_div_u_div_v!(A::Matrix{Float64}, coord::Matrix{Float64},
	                  a::Float64)
    b, area = barycentric(coord)
    c = 4 * area * a
    for j = 1:3
        j₋ = PRED[j]
        for i = 1:j
            i₋ = PRED[i]
            A[i,j]     = c * b[1,j₋] * b[1,i₋]
            A[i+3,j+3] = c * b[2,j₋] * b[2,i₋]
        end
    end
    for j = 1:2
        for i = j+1:3
            A[i,j]     = A[j,i] 
            A[i+3,j+3] = A[j+3,i+3]
        end
    end
    for j = 1:3
        j₋ = PRED[j]
        for i = 1:3
            i₋ = PRED[i]
            A[i+3,j] = c * b[1,j₋] * b[2,i₋]
        end
    end
    for j = 1:3
        for i = 1:3
            A[i,j+3] = A[j+3,i]
        end
    end
end

function ∫∫a_div_u_div_v!(A::Matrix{Float64}, coord::Matrix{Float64},
	                  a::Function, params...)
    Q = 0.0
    for j = 1:3
	p = j + 3
	x = coord[1,p] # coordinates of jth midpoint
	y = coord[2,p] 
	Q += a(x, y, params...)
    end
    Q /= 3
    ∫∫a_div_u_div_v!(A, coord, Q)
end

∫∫a_div_u_div_v!() = ("Triangle 6", [4, 5, 6])

function ∫∫μ_∇u_colon_∇v!(B::Matrix{Float64}, coord::Matrix{Float64},
	                  μ::Float64)
    b, area = barycentric(coord)
    a = 4 * μ * area
    for j = 1:3
	j₋ = PRED[j]
	for i = j:3
	    i₋ = PRED[i]
	    B[i,j] = a * ( b[1,j₋] * b[1,i₋] + b[2,j₋] * b[2,i₋] )
	    B[i+3,j+3] = B[i,j]
	end
    end
    for j = 2:3
	for i = 1:j-1
	    B[i,j] = B[j,i]
	    B[i+3,j+3] = B[j+3,i+3]
	end
    end
    B[4:6,1:3] .= 0.0
    B[1:3,4:6] .= 0.0
end

function ∫∫μ_∇u_colon_∇v!(C::Matrix{Float64}, coord::Matrix{Float64},
	                  μ::Function, params...)
    Q = 0.0
    for j = 1:3
	x = coord[1,j+3] # coordinates of jth midpoint
	y = coord[2,j+3] 
	Q += μ(x, y, params...)
    end
    Q /= 3
    ∫∫μ_∇u_colon_∇v!(C, coord, Q)
end

∫∫μ_∇u_colon_∇v!() = ("Triangle 6", [4, 5, 6])
    
function correction!(D::Matrix{Float64}, coord::Matrix{Float64},
	                ∇μ::SVector)  
    b, area = barycentric(coord)
    c = (2/3) * area
    for j = 1:3
	for i = 1:3
	    i₋ = PRED[i]
	    D[i,j] = 0.0
	    D[i,j+3] = c * ( ∇μ[1] * b[2,i₋] - ∇μ[2] * b[1,i₋] )
	    D[j+3,i] = D[i,j+3]
	    D[i+3,j+3] = 0.0
	end
    end
end

function correction!(D::Matrix{Float64}, coord::Matrix{Float64},
	             ∇μ::Function, params...)  
    b, area = barycentric(coord)
    c = (2/3) * area
    for j = 1:3
	x = coord[1,j+3]
	y = coord[2,j+3]
	v = ∇μ(x, y, params...) # value at jth midpoint
	for i = 1:3
	    i₋ = PRED[i]
	    D[i,j] = 0.0
	    D[i,j+3] = c * ( v[1] * b[2,i₋] - v[2] * b[1,i₋] )
	    D[j+3,i] = D[i,j+3]
	    D[i+3,j+3] = 0.0
	end
    end
end

correction!() = ("Triangle 6", [4, 5, 6])

function ∫∫f_dot_v!(v::Vector{Float64}, coord::Matrix{Float64},
                    f::SVector{2,Float64})
    b, area = barycentric(coord)
    v[1:3] .= area * f[1] / 3
    v[4:6] .= area * f[2] / 3
end

function ∫∫f_dot_v!(v::Vector{Float64}, coord::Matrix{Float64},
                    f::Function, params...)
    b, area = barycentric(coord)
    c = area / 3
    for i = 1:3
	p = i + 3
	x = coord[1,p]
	y = coord[2,p]
	fval = f(x, y, params...)
        v[i] = c * fval[1] 
        v[p] = c * fval[2]
    end
end

∫∫f_dot_v!() = ("Triangle 6", [4, 5, 6])

function ∫g_dot_v!(loadvec::Vector{Float64}, coord::Matrix{Float64},
                   g::SVector{2,Float64})
    x1, y1 = coord[1,1], coord[2,1]
    x2, y2 = coord[1,2], coord[2,2]
    lenE = hypot(x1-x2, y1-y2)
    loadvec[1] = lenE * g[1]
    loadvec[2] = lenE * g[2]
end

function ∫g_dot_v!(loadvec::Vector{Float64}, coord::Matrix{Float64},
	g::Function, params...)
    x1, y1 = coord[1,1], coord[2,1]
    x2, y2 = coord[1,2], coord[2,2]
    lenE = hypot(x1-x2, y1-y2)
    x3, y3 = coord[1,3], coord[2,3]
    loadvec .= lenE * g(x3, y3, params...)
end

∫g_dot_v!() = ("Line 3", [3])

function error_norms(u1h::Vector{Float64}, u2h::Vector{Float64}, u::Function,
                 ∇u::Function, dof::DegreesOfFreedom, level::Int64, params...)
    mesh = dof.mesh
    L2err, H1err = 0.0, 0.0
    ψ = Vector{Float64}(undef, 3)
    b = Matrix{Float64}(undef, 2, 3)
    cntrd = Vector{Float64}(undef, 2)
    ∇uh = Matrix{Float64}(undef, 2, 2)
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
            fill!(∇uh, 0.0)
            for j = eachindex(elt_global_dof)
                r = elt_global_dof[j]
		j₋ = PRED[j]
		∇uh[1,1] += u1h[r] * (-2 * b[1,j₋])
		∇uh[2,1] += u1h[r] * (-2 * b[2,j₋])
		∇uh[1,2] += u2h[r] * (-2 * b[1,j₋])
		∇uh[2,2] += u2h[r] * (-2 * b[2,j₋])
            end
            Σ_L2, Σ_H1 = 0.0, 0.0
            for k in eachindex(x)
                nc_shape_func!(ψ, x[k], y[k], b, cntrd)
                u1h_k = 0.0
                u2h_k = 0.0
                for j in eachindex(elt_global_dof)
                    r = elt_global_dof[j]
                    u1h_k += u1h[r] * ψ[j]
                    u2h_k += u2h[r] * ψ[j]
                end
                u_k = u(x[k], y[k], params...)
                ∇u_k = ∇u(x[k], y[k], params...)
                Σ_L2 += ( u1h_k - u_k[1] )^2 + ( u2h_k - u_k[2] )^2
                for p = 1:4
                    Σ_H1 += ( ∇uh[p] - ∇u_k[p] )^2
                end
            end
            L2err += (area / lastindex(x)) * Σ_L2
            H1err += (area / lastindex(x)) * Σ_H1
        end
    end
    L2err = sqrt(L2err)
    H1err = sqrt(H1err)
    return L2err, H1err
end

function L2error(u1h::Vector{Float64}, u2h::Vector{Float64}, u::Function,
                 dof::DegreesOfFreedom, level::Int64, params...)
    mesh = dof.mesh
    L2err = 0.0
    ψ = Vector{Float64}(undef, 3)
    b = Matrix{Float64}(undef, 2, 3)
    cntrd = Vector{Float64}(undef, 2)
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
            Σ_L2 = 0.0
            for k in eachindex(x)
                nc_shape_func!(ψ, x[k], y[k], b, cntrd)
                u1h_k = 0.0
                u2h_k = 0.0
                for j in eachindex(elt_global_dof)
                    r = elt_global_dof[j]
                    u1h_k += u1h[r] * ψ[j]
                    u2h_k += u2h[r] * ψ[j]
                end
                u_k = u(x[k], y[k], params...)
                Σ_L2 += ( u1h_k - u_k[1] )^2 + ( u2h_k - u_k[2] )^2
            end
            L2err += (area / lastindex(x)) * Σ_L2
        end
    end
    L2err = sqrt(L2err)
    return L2err
end
end # module
