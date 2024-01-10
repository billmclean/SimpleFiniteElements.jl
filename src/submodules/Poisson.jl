module Poisson

import ..SUCC, ..FEMesh, ..DegreesOfFreedom
import ..FEM: assemble_matrix, assemble_vector, prepare_assembly
import ..Utils: barycentric, triangle_quad, centroid, 
                quadrature_points, shape_func!, barycentric!
import ..MeshGen: elt_node_coord!
import LinearAlgebra: cholesky, lmul!
import StaticArrays: SA

function poisson_soln(dof::DegreesOfFreedom,
	bilinear_forms::Dict, linear_funcs::Dict)
    mesh = dof.mesh
    A_free, A_fix = assemble_matrix(dof, bilinear_forms)
    b_free, u_fix = assemble_vector(dof, linear_funcs)
    b = b_free - A_fix * u_fix
    CF = cholesky(A_free)
    u_free = CF \ b
    uh = [ u_free; u_fix ]
    return uh
end

"""
    ∫∫a_∇u_dot_∇v_dA!(A, coord, dof)

Computes the 3×3 element stiffness matrix `A` for a triangle with
vertices `coord[:,j]`, where `a(x,y) = coef`.
"""
function ∫∫a_∇u_dot_∇v!(A::Matrix{Float64}, coord::Matrix{Float64},
                      coef::Float64) 
    b, area = barycentric(coord)
    for q=1:3
        for p = q:3
            A[p,q] = area * ( b[1,p] * b[1,q] + b[2,p] * b[2,q] )
        end
    end
    for q = 2:3
        for p = 1:q-1
            A[p,q] = A[q,p]
        end
    end
    lmul!(coef, A)
end

function ∫∫a_∇u_dot_∇v!(A::Matrix{Float64}, coord::Matrix{Float64},
	              coef::Function, params...)
    avg = triangle_quad(coord, coef, params...)
    ∫∫a_∇u_dot_∇v!(A, coord, avg)
end

∫∫a_∇u_dot_∇v!() = ("Triangle 3", [1, 2, 3])

function ∫∫c_u_v!(M::Matrix{Float64}, coord::Matrix{Float64}, coef::Float64)
    b, area = barycentric(coord)
    fill!(M, area/12)
    for j = 1:3
	M[j,j] *= 2
    end
    lmul!(coef, M)
end

function ∫∫c_u_v!(M::Matrix{Float64}, coord::Matrix{Float64}, 
	coef::Function, params...)
    b, area = barycentric(coord)
    x1, y1 = coord[1,1], coord[2,1]
    x2, y2 = coord[1,2], coord[2,2]
    x3, y3 = coord[1,3], coord[2,3]
    c1 = coef(x1, y1, params...)
    c2 = coef(x2, y2, params...)
    c3 = coef(x3, y3, params...)
    D = area / 60
    M[1,1] = D * ( 6c1 + 2c2 + 2c3 )
    M[2,1] = D * ( 2c1 + 2c2 +  c3 )
    M[3,1] = D * ( 2c1 +  c2 + 2c3 )
    M[1,2] = M[2,1]
    M[2,2] = D * ( 2c1 + 6c2 + 2c3 )
    M[3,2] = D * (  c1 + 2c2 + 2c3 )
    M[1,3] = M[3,1]
    M[2,3] = M[3,2]
    M[3,3] = D * ( 2c1 + 2c2 + 6c3 )
end

∫∫c_u_v!() = ("Triangle 3", [1, 2, 3])

function ∫c_u_v!(M::Matrix{Float64}, coord::Matrix{Float64}, coef::Float64)
    h = hypot(coord[1,2]-coord[1,1], coord[2,2]-coord[2,1])
    M[1,1] = h / 3
    M[2,1] = h / 6
    M[1,2] = M[2,1]
    M[2,2] = M[1,1]
    lmul!(coef, M)
end

function ∫c_u_v!(M::Matrix{Float64}, coord::Matrix{Float64}, 
	coef::Function, params...)
    x1, y1 = coord[1,1], coord[2,1]
    x2, y2 = coord[1,2], coord[2,2]
    h = hypot(x1-x2, y1-y2)
    c1 = coef(x1, y1, params...)
    c2 = coef(x2, y2, params...)
    D = h / 12
    M[1,1] = D * ( 3c1 +  c2 )
    M[2,1] = D * (  c1 +  c2 )
    M[1,2] = M[2,1]
    M[2,2] = D * (  c1 + 3c2 )
end

∫c_u_v!() = ("Line 2", [1, 2])

function ∫∫f_v!(v::Vector{Float64}, coord::Matrix{Float64}, f::Float64)
    b, area = barycentric(coord)
    fill!(v, f * area / 3)
end

function ∫∫f_v!(v::Vector{Float64}, coord::Matrix{Float64}, 
	f::Function, params...)
    b, area = barycentric(coord)
    f1 = f(coord[1,1], coord[2,1], params...)
    f2 = f(coord[1,2], coord[2,2], params...)
    f3 = f(coord[1,3], coord[2,3], params...)
    c = area/12
    v[1] = c * ( 2*f1 +   f2 +   f3 )
    v[2] = c * (   f1 + 2*f2 +   f3 )
    v[3] = c * (   f1 +   f2 + 2*f3 )
end

∫∫f_v!() = ("Triangle 3", [1, 2, 3])

function ∫g_v!(v::Vector{Float64}, coord::Matrix{Float64}, g::Float64)
    lenE = hypot(coord[1,2]-coord[1,1], coord[2,2]-coord[2,1])
    fill!(v, g*lenE/2)
end

function ∫g_v!(v::Vector{Float64}, coord::Matrix{Float64}, 
	g::Function, params...)
    x1, y1 = coord[1,1], coord[2,1]
    g1 = g(x1, y1, params...)
    x2, y2 = coord[1,2], coord[2,2]
    g2 = g(x2, y2, params...)
    lenE = hypot(x2-x1, y2-y1)
    v[1] = (lenE/6) * ( 2g1 +  g2 )
    v[2] = (lenE/6) * (  g1 + 2g2 )
end

∫g_v!() = ("Line 2", [1, 2])

function error_norms(uh::Vector{Float64}, u::Function, ∇u::Function, 
	             dof::DegreesOfFreedom, level::Int64, params...)
    mesh = dof.mesh
    L2err, H1err = 0.0, 0.0
    ψ = Vector{Float64}(undef, 3)
    b = Matrix{Float64}(undef, 2, 3)
    cntrd = Vector{Float64}(undef, 2)
    ∇uh = similar(cntrd)
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
            fill!(∇uh, 0.0) # ∇uh is constant on each element
	    for j = eachindex(elt_global_dof)
                r = elt_global_dof[j]
		∇uh[1] += uh[r] * b[1,j]
		∇uh[2] += uh[r] * b[2,j]
	    end
            Σ_L2, Σ_H1 = 0.0, 0.0
            for k in eachindex(x)
                shape_func!(ψ, x[k], y[k], b, cntrd)
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
