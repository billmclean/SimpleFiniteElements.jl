module Utils

import ..DegreesOfFreedom, ..FEMesh, ..PRED, ..SUCC
import Base
import LinearAlgebra
import LinearAlgebra: cholesky, mul!, ldiv!
import SparseArrays: SparseMatrixCSC
import SuiteSparse

"""
    b, area = barycentric(coord)

For a planar triangle with vertices `z_j=coord[:,j]`
compute vectors `bⱼ= b[:,j]` such that if

    x = λ₁ z₁ + λ₂ z₂ + λ₃ z₃

then the barycentric coordinates of `x` are given by

    λₖ(x) = 1/3 + bₖ ⋅ ( x - centroid )

Also compute the area of the triangle.
"""
function barycentric(coord::Matrix{Float64})
    b = Matrix{Float64}(undef, 2, 3)
    b[1,1] = coord[1,1] - coord[1,3]
    b[2,1] = coord[2,1] - coord[2,3]
    b[1,2] = coord[1,2] - coord[1,3]
    b[2,2] = coord[2,2] - coord[2,3]
    # overwrite B[1:2,1:2] with its inverse transpose
    area = b[1,1] * b[2,2] - b[1,2] * b[2,1]
    b[1,1], b[2,2] =  b[2,2]/area,  b[1,1]/area
    b[1,2], b[2,1] = -b[2,1]/area, -b[1,2]/area

    b[1,3] = -b[1,1] - b[1,2]
    b[2,3] = -b[2,1] - b[2,2]
    area = abs(area) / 2
    return b, area
end

function barycentric!(b::Matrix{Float64}, cntrd::Vector{Float64}, 
	coord::Matrix{Float64})
    b[1,1] = coord[1,1] - coord[1,3]
    b[2,1] = coord[2,1] - coord[2,3]
    b[1,2] = coord[1,2] - coord[1,3]
    b[2,2] = coord[2,2] - coord[2,3]
    # overwrite B[1:2,1:2] with its inverse transpose
    area = b[1,1] * b[2,2] - b[1,2] * b[2,1]
    b[1,1], b[2,2] =  b[2,2]/area,  b[1,1]/area
    b[1,2], b[2,1] = -b[2,1]/area, -b[1,2]/area

    b[1,3] = -b[1,1] - b[1,2]
    b[2,3] = -b[2,1] - b[2,2]
    area = abs(area) / 2
    cntrd[1] = ( coord[1,1] + coord[1,2] + coord[1,3] ) / 3
    cntrd[2] = ( coord[2,1] + coord[2,2] + coord[2,3] ) / 3
    return area
end

function centroid(coord::Matrix)
    cntrd = sum(coord, dims=2) / 3 
    return cntrd[:]
end

function shape_func!(ψ::AbstractVector{Float64}, x::Float64, y::Float64, 
	b::Matrix{Float64}, cntrd::Vector{Float64})
    Δx = x - cntrd[1]
    Δy = y - cntrd[2]
    for j = 1:3
	ψ[j] = 1/3 + b[1,j] * Δx + b[2,j] * Δy
    end
end

function shape_func(x, y, b, cntrd) 
    ψ = Vector{Float64}(undef, 3)
    shape_func!(ψ, x, y, b, cntrd)
    return ψ
end

function nc_shape_func!(ψ::AbstractVector{Float64}, x::Float64, y::Float64, 
	b::Matrix{Float64}, cntrd::Vector{Float64})
    shape_func!(ψ, x, y, b, cntrd)
    temp = ψ[3]
    ψ[3] = 1 - 2ψ[2]
    ψ[2] = 1 - 2ψ[1]
    ψ[1] = 1 - 2temp
end

function nc_shape_func(x, y, b, cntrd) 
    ψ = Vector{Float64}(undef, 3)
    nc_shape_func!(ψ, x, y, b, cntrd)
    return ψ
end

function diameter(coord::Matrix{Float64})
    d = 0.0
    for j = 1:3
        pj = PRED[j]
        sj = SUCC[j]
        Δx = coord[1,pj] - coord[1,sj]
        Δy = coord[2,pj] - coord[2,sj]
        side_len = hypot(Δx, Δy)
        d = max(side_len, d)
    end
    return d
end

function gmsh2pyplot(dof::DegreesOfFreedom)
    mesh = dof.mesh
    node_tag = dof.node_tag
    x = Vector{Float64}(undef, lastindex(node_tag))
    y = similar(x)
#    for nd = 1:num_mesh_nodes
    for nd in eachindex(node_tag)
	j = dof.node_tag_to_dof[nd]
	x[j] = mesh.coord[1,nd]
	y[j] = mesh.coord[2,nd]
    end
    physical_groups = mesh.gmodel.physical_groups
    num_triangles = 0
    for (name, dim_tag) in physical_groups
	dim, tag = dim_tag
	if dim == 2
	    et = mesh.elt_type_in[name]
	    @assert mesh.elt_properties[et].element_name == "Triangle 3"
	    num_triangles += length(mesh.elt_tags_in[name])
	end
    end
    triangles = Matrix{Int64}(undef, (num_triangles,3))
    i = 0
    for (name, elt_node_tags) in mesh.elt_node_tags_in
	dim, tag = physical_groups[name]
	if dim == 2
	    elt_tags = mesh.elt_tags_in[name]
	    m = length(elt_tags)
	    elt_node_tags = reshape(elt_node_tags, (3,m))
	    for k = 1:m
		i += 1
		for j = 1:3
		    nd = elt_node_tags[j,k]
		    l = dof.node_tag_to_dof[nd]
  		    triangles[i,j] = l - 1
		end
	    end
	end
    end
    return x, y, triangles
end

"""
    triangle_quad(coord, f, params...)

Returns a quadrature approximation to ∫∫f / area, where the integration
domain is the triangle with vertices `coord[:,j]` for `j = 1, 2, 3`.  This
equal-weight rule evaluates `f(x, y, params...)` at the midpoints of the 
triangle edges, and integrates quadratic polynomials exactly.
"""
function triangle_quad(coord::Matrix{Float64}, f::Function, params...)
    Q = 0.0
    for j = 1:3
        sj = SUCC[j]
        x = ( coord[1,j] + coord[1,sj] ) / 2
        y = ( coord[2,j] + coord[2,sj] ) / 2
        Q += f(x, y, params...)
    end
    Q /= 3
    return Q
end

"""
    quadrature_points(level, coord) -> x, y

Generates an equal-weight, composite quadrature rule by subdividing the triangle
with vertices `coord[:,k]` into `4^level` congruent sub-triangles, each 
containing `3` quadrature points.  The basic `3`-point rule is exact for
quadratic polynomials so the composite rule is accurate to order `h^3`
for `h = 2^(-level)`.  The weights are all equal to `1/M` where 
`M = 3 × 4^level` is the number of quadrature points.
"""
function quadrature_points(level::Integer, coord::Matrix{T}
                          ) where T <: AbstractFloat
    N = 4^level
    n = 2^level
    x = Matrix{T}(undef, 3*N, 2)
    fill!(x, NaN)
    k = 0
    Δ₁ = ( coord[:,2] - coord[:,1] ) / n
    Δ₂ = ( coord[:,3] - coord[:,1] ) / n
    l = 0
    one_sixth = one(T) / 6
    a₁ = Vector{T}(undef, 2)
    a₂ = similar(a₁)
    a₃ = similar(a₁)
    δ₁ = similar(a₁)
    δ₂ = similar(a₁)
    for i = 1:n
        δ₁ .= (i-1) * Δ₁
        for j = 1:n - i + 1
            δ₂ .= (j-1) * Δ₂
            a₁ .= coord[:,1] .+ δ₁ .+ δ₂
            a₂ .= a₁ .+ Δ₁
            a₃ .= a₁ .+ Δ₂
            x[l+1,:] .= one_sixth * ( 4a₁ .+  a₂ .+  a₃ )
            x[l+2,:] .= one_sixth * (  a₁ .+ 4a₂ .+  a₃ )
            x[l+3,:] .= one_sixth * (  a₁ .+  a₂ .+ 4a₃ )
            l += 3
        end
    end
    for i = 1:n-1
        δ₁ .= i * Δ₁
        for j = 1:n - i
            δ₂ .= j * Δ₂
            a₁ .= coord[:,1] .+ δ₁ .+ δ₂
            a₂ .= a₁ .- Δ₁
            a₃ .= a₁ .- Δ₂
            x[l+1,:] .= one_sixth * ( 4a₁ .+  a₂ .+  a₃ )
            x[l+2,:] .= one_sixth * (  a₁ .+ 4a₂ .+  a₃ )
            x[l+3,:] .= one_sixth * (  a₁ .+  a₂ .+ 4a₃ )
            l += 3
        end
    end
    return x[:,1], x[:,2]
end

"""
    SchurComplement

The linear operator `x → Sx` where

    S = A₁₁ - A₁₂ inv(A₂₂) A₂₁

"""
struct SchurComplement 
    A11 :: SparseMatrixCSC{Float64}
    A12 :: SparseMatrixCSC{Float64}
    CA  :: SuiteSparse.CHOLMOD.Factor{Float64}
    z1  :: Vector{Float64}
    z2  :: Vector{Float64}
end

function SchurComplement(A :: SparseMatrixCSC, n1::Int64)
    n = size(A, 1)
    size(A, 2) == n || throw(DimensionMismatch("A must be a square matrix"))
    n2 = n - n1
    (0 < n2 < n) || error("We require 1 ≤ n1 ≤ size(A, 1)")
    A11 = A[1:n1,1:n1]
    A12 = A[1:n1,n1+1:n]
    A22 = A[n1+1:n,n1+1:n]
    CA = cholesky(A22)
    z1 = Vector{Float64}(undef, n2)
    z2 = Vector{Float64}(undef, n2)
    return SchurComplement(A11, A12, CA, z1, z2)
end

Base.size(S::SchurComplement) = size(S.A11)
Base.size(S::SchurComplement, i::Int64) = size(S.A11, i)
Base.eltype(S::SchurComplement) = Float64

function LinearAlgebra.mul!(y::AbstractVector{Float64}, S::SchurComplement, 
              x::Vector{Float64})
    A11, A12, CA, z1, z2 = S.A11, S.A12, S.CA, S.z1, S.z2
    n1 = size(A11, 1)
    length(y) == n1 || throw(DimensionMismatch)
    length(x) == n1 || throw(DimensionMismatch)
    A21 = transpose(A12)
    mul!(z1, A21, x)           # z1 = A_21 x
    z2 .= CA \ z1              # z2 = A_22 \ ( A_21 x )
    mul!(y, A12, z2)           # y = A_12 inv(A_22) A_21 x
    mul!(y, A11, x, 1.0, -1.0) # y = A_11 x - y
end

function LinearAlgebra.mul!(Y::AbstractMatrix{Float64}, S::SchurComplement,
	X::AbstractMatrix{Float64})
    A11, A12, CA = S.A11, S.A12, S.CA
    n1, n2 = size(S)
    size(X, 1) == n1 || throw(DimensionMismatch)
    size(Y, 1) == n1 || throw(DimensionMismatch)
    n3 = size(X, 2)
    size(Y, 2) == n3 || throw(DimensionMismatch)
    A21 = transpose(A12)
    Z1 = A21 * X
    Z2 = CA \ Z1
    mul!(Y, A12, Z2)
    mul!(Y, A11, X, 1.0, -1.0)
end

function Base.Matrix(S::SchurComplement)
    n1, n2 = size(S)
    n1 == n2 || throw(DimensionMismatch)
    e = zeros(n1)
    Smat = zeros(n1, n1)
    e[1] = 1.0
    col1 = view(Smat, :, 1)
    mul!(col1, S, e)
    for j = 2:n1
	e[j-1] = 0.0
	e[j]   = 1.0
	colj = view(Smat, :, j)
	mul!(colj, S, e)
    end
    return Smat
end

end # module
