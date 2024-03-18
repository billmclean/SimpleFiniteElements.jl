module Elasticity

import ..FEMesh, ..DegreesOfFreedom
import ..FEM: assemble_vector, assemble_matrix, prepare_assembly
import StaticArrays: SVector, SA
import ..Utils: barycentric, gmsh2pyplot, triangle_quad, centroid,
                quadrature_points, shape_func!, barycentric!
import ..MeshGen: elt_node_coord!
import PyPlot: figure, quiver, axis, xlabel, ylabel, plot
import Statistics: mean
import LinearAlgebra: cholesky

const ELT_DOF = Dict(Int32(1) => [1, 2],
                     Int32(2) => [1, 2, 3])

function elasticity_soln(dof::DegreesOfFreedom, 
	bilinear_forms::Dict, linear_funcs::Dict)
    mesh = dof.mesh
    b_free, u_fix = assemble_vector(dof, linear_funcs, 2)
    A_free, A_fix = assemble_matrix(dof, bilinear_forms, 2)
    b = b_free - A_fix * u_fix
    F = cholesky(A_free)
    u_free = F \ b
    num_free, num_fixed = dof.num_free, dof.num_fixed
    u1h = [ u_free[1:num_free]; u_fix[1:num_fixed] ]
    u2h = [ u_free[num_free+1:2*num_free]; u_fix[num_fixed+1:2*num_fixed] ]
    return u1h, u2h
end

function elasticity_soln(dof::DegreesOfFreedom, 
	bilinear_forms::Vector{<:Tuple},
	linear_funcs::Vector{<:Tuple})
    mesh = dof.mesh
    b_free, u_fix = assemble_vector(dof, linear_funcs, 2)
    A_free, A_fix = assemble_matrix(dof, bilinear_forms, 2)
    b = b_free - A_fix * u_fix
    u_free = A_free \ b
    num_free, num_fixed = dof.num_free, dof.num_fixed
    u1h = [ u_free[1:num_free]; u_fix[1:num_fixed] ]
    u2h = [ u_free[num_free+1:2*num_free]; u_fix[num_fixed+1:2*num_fixed] ]
    return u1h, u2h
end

function visualise_soln(dof::DegreesOfFreedom,
	u1h::Vector{Float64}, u2h::Vector{Float64}, scale::Float64, figno=1)
    mesh = dof.mesh
    x, y, triangles = gmsh2pyplot(dof)
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
    function add_point!(x, y, x_ref, y_ref, nd)
	    j = node_tag_to_dof[nd]
	    push!(x, coord[1,nd]+u1h[j]/scale)
	    push!(y, coord[2,nd]+u2h[j]/scale)
	    push!(x_ref, coord[1,nd])
	    push!(y_ref, coord[2,nd])
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
	    nd1 = elt_node_tags[1,1]
	    add_point!(x, y, x_ref, y_ref, nd1)
	    nd2 = elt_node_tags[2,1]
	    add_point!(x, y, x_ref, y_ref, nd2)
	    for j = 2:num_elts
		nd1 = elt_node_tags[1,j]
		if nd1 ≠ elt_node_tags[j-1]
		    push!(x, NaN)
		    push!(y, NaN)
		    push!(x_ref, NaN)
		    push!(y_ref, NaN)
		    add_point!(x, y, x_ref, y_ref, nd1)
		end
		nd2 = elt_node_tags[2,j]
		add_point!(x, y, x_ref, y_ref, nd2)
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

function ∫∫λ_div_u_div_v!(A::Matrix{Float64}, coord::Matrix{Float64}, 
	                  λ::Float64)
    b, area = barycentric(coord)
    c = area * λ
    for j = 1:3
	for i = 1:j
	    A[i,j]     = c * b[1,j] * b[1,i]
	    A[i+3,j+3] = c * b[2,j] * b[2,i]
	end
    end
    for j = 1:2
	for i = j+1:3
	    A[i,j]     = A[j,i]
	    A[i+3,j+3] = A[j+3,i+3]
	end
    end
    for j = 1:3
	for i = 1:3
	    A[i+3,j] = c * b[1,j] * b[2,i]
	end
    end
    for j = 1:3
	for i = 1:3
	    A[i,j+3] = A[j+3,i]
	end
    end
end

function ∫∫λ_div_u_div_v!(A::Matrix{Float64}, coord::Matrix{Float64}, 
	                  λ::Function, params...)
    Q = triangle_quad(coord, λ, params...)
    ∫∫λ_div_u_div_v!(A, coord, Q)
end

∫∫λ_div_u_div_v!() = ("Triangle 3", [1, 2, 3])

function ∫∫2μ_εu_εv!(B::Matrix{Float64}, coord::Matrix{Float64},
	            μ::Float64)
    b, area = barycentric(coord)
    c = area * μ
    for j = 1:3
	for i = 1:j
	    B[i,j]     = c * ( 2*b[1,j]*b[1,i] +   b[2,j]*b[2,i] )
	    B[i+3,j+3] = c * (   b[1,j]*b[1,i] + 2*b[2,j]*b[2,i] )
	end
    end
    for j = 1:2
	for i = j+1:3
	    B[i,j] = B[j,i]
	    B[i+3,j+3] = B[j+3,i+3]
	end
    end
    for j = 1:3
	for i = 1:3
	    B[i+3,j] = c * b[2,j] * b[1,i]
	end
    end
    for j = 1:3
	for i = 1:3
	    B[i,j+3] = B[j+3,i]
	end
    end
end

function ∫∫2μ_εu_εv!(B::Matrix{Float64}, coord::Matrix{Float64},
	            μ::Function, params...)
    Q = triangle_quad(coord, μ, params...)
    ∫∫2μ_εu_εv!(B, coord, Q)
end

∫∫2μ_εu_εv!() = ("Triangle 3", [1, 2, 3])

function ∫∫c_u_dot_v!(M::Matrix{Float64}, coord::Matrix{Float64},
	              c::Float64)
    b, area = barycentric(coord)
    M[1:3,1:3] .= area / 12
    for j = 1:3
	M[j,j] *= 2
    end
    M[1:3,4:6] .= 0.0
    M[4:6,1:3] .= 0.0
    M[4:6,4:6] .= M[1:3,1:3]
end

∫∫c_u_dot_v!() = ("Triangle 3", [1, 2, 3])

function ∫∫f_dot_v!(v::Vector{Float64}, coord::Matrix{Float64}, 
	            f::SVector{2,Float64})
    b, area = barycentric(coord)
    v[1:3] .= area * f[1] / 3
    v[4:6] .= area * f[2] / 3
end

function ∫∫f_dot_v!(v::Vector{Float64}, coord::Matrix{Float64}, 
                    f::Function, params...)
    b, area = barycentric(coord)
    f1 = f(coord[1,1], coord[2,1], params...)
    f2 = f(coord[1,2], coord[2,2], params...)
    f3 = f(coord[1,3], coord[2,3], params...)
    c = area/12
    v[1] = c * ( 2*f1[1] +   f2[1] +   f3[1] )
    v[2] = c * (   f1[1] + 2*f2[1] +   f3[1] )
    v[3] = c * (   f1[1] +   f2[1] + 2*f3[1] )
    v[4] = c * ( 2*f1[2] +   f2[2] +   f3[2] )
    v[5] = c * (   f1[2] + 2*f2[2] +   f3[2] )
    v[6] = c * (   f1[2] +   f2[2] + 2*f3[2] )
end

∫∫f_dot_v!() = ("Triangle 3", [1, 2, 3])

function ∫g_dot_v!(v::Vector{Float64}, coord::Matrix{Float64},
	           g::SVector{2,Float64})
    x1, y1 = coord[1,1], coord[2,1]
    x2, y2 = coord[1,2], coord[2,2]
    lenE = hypot(x1-x2, y1-y2)
    C = lenE / 2
    v[1] = v[2] = C * g[1]
    v[3] = v[4] = C * g[2]
end

function ∫g_dot_v!(v::Vector{Float64}, coord::Matrix{Float64}, 
	           g::Function, params...)
    x1, y1 = coord[1,1], coord[2,1]
    g1 = g(x1, y1, params...)
    x2, y2 = coord[1,2], coord[2,2]
    g2 = g(x2, y2, params...)
    lenE = hypot(x1-x2, y1-y2)
    C = lenE / 6
    v[1] = C * ( 2g1[1] +  g2[1] )
    v[2] = C * (  g1[1] + 2g2[1] )
    v[3] = C * ( 2g1[2] +  g2[2] )
    v[4] = C * (  g1[2] + 2g2[2] )
end

∫g_dot_v!() = ("Line 2", [1, 2])

function fundamental_soln(x::SVector{2,Float64}, y::SVector{2,Float64},
	λ::Float64, μ::Float64)
    r = x - y
    norm_r = hypot(r[1], r[2])
    C1 = (3μ+λ) / ( 4π * μ * (2μ+λ) )
    C2 =  (μ+λ) / ( 4π * μ * (2μ+λ) )
    lg = log(1/norm_r)
    ω = r / norm_r
    G11 = C1 * lg + C2 * ω[1] * ω[1]
    G21 = G12 = C2 * ω[1] * ω[2]
    G22 = C1 * lg + C2 * ω[2] * ω[2]
    return SA[ G11 G12
	       G21 G22 ]
end

function L2error(u1h::Vector{Float64}, u2h::Vector{Float64}, u::Function, 
                 dof::DegreesOfFreedom, level::Int64, params...)
    mesh = dof.mesh
    L2err = 0.0
    ψ = Vector{Float64}(undef, 3)
    b = Matrix{Float64}(undef, 2, 3)
    cntrd = Vector{Float64}(undef, 2)
    for (name, elt_type) in mesh.elt_type_in
        # Ignore unless element type is "Triangle 3".
        if elt_type ≠ 2
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
                shape_func!(ψ, x[k], y[k], b, cntrd)
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
function error_norms(u1h::Vector{Float64}, u2h::Vector{Float64}, u::Function, 
                 ∇u::Function, dof::DegreesOfFreedom, level::Int64, params...)
    mesh = dof.mesh
    L2err, H1err = 0.0, 0.0
    ψ = Vector{Float64}(undef, 3)
    b = Matrix{Float64}(undef, 2, 3)
    cntrd = Vector{Float64}(undef, 2)
    ∇uh = Matrix{Float64}(undef, 2, 2)
    for (name, elt_type) in mesh.elt_type_in
        # Ignore unless element type is "Triangle 3".
        if elt_type ≠ 2
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
		∇uh[1,1] += u1h[r] * b[1,j]
		∇uh[2,1] += u1h[r] * b[2,j]
		∇uh[1,2] += u2h[r] * b[1,j]
		∇uh[2,2] += u2h[r] * b[2,j]
	    end
            Σ_L2, Σ_H1 = 0.0, 0.0
            for k in eachindex(x)
                shape_func!(ψ, x[k], y[k], b, cntrd)
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

end # module
