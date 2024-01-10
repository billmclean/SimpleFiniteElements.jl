import SimpleFiniteElements: PRED, SUCC
import SimpleFiniteElements.Elasticity: ∫∫λ_div_u_div_v!, ∫∫2μ_εu_εv!, 
                                        ∫∫f_dot_v!
import SimpleFiniteElements.NonConformingElasticity as NCE
import SimpleFiniteElements.Utils: barycentric, centroid, shape_func, 
				   quadrature_points, nc_shape_func
import StaticArrays: SA
import LinearAlgebra: dot, issymmetric

function vector_shape_func(x::Float64, y::Float64, b::Matrix{Float64},
	cntrd::Vector{Float64})
    ψ = zeros(2, 6)
    ψ[1,1:3] = shape_func(x, y, b, cntrd)
    ψ[2,4:6] .= ψ[1,1:3]
    return ψ
end

function nc_vector_shape_func(x::Float64, y::Float64, b::Matrix{Float64},
	cntrd::Vector{Float64})
    ψ = zeros(2, 6)
    ψ[1,1:3] = nc_shape_func(x, y, b, cntrd)
    ψ[2,4:6] = ψ[1,1:3]
    return ψ
end

function shape_strain!(εψ::Array{Float64}, x::Float64, y::Float64, 
	               b::Matrix{Float64}, cntrd::Vector{Float64}, h::Float64)
    ψ₁₊ = vector_shape_func(x+h, y, b, cntrd)
    ψ₁₋ = vector_shape_func(x-h, y, b, cntrd)
    ψ₂₊ = vector_shape_func(x, y+h, b, cntrd)
    ψ₂₋ = vector_shape_func(x, y-h, b, cntrd)
    for j = 1:6
	εψ[1,:,j] = ( ψ₁₊[:,j] - ψ₁₋[:,j] ) / (2h)
	εψ[2,:,j] = ( ψ₂₊[:,j] - ψ₂₋[:,j] ) / (2h)
	εψ[1,2,j] = ( εψ[1,2,j] + εψ[2,1,j] ) / 2
	εψ[2,1,j] = εψ[1,2,j]
    end
end

function nc_shape_gradient!(∇ψ::Array{Float64}, x::Float64, y::Float64,
	                 b::Matrix{Float64}, cntrd::Vector{Float64}, h::Float64)
    ψ₁₊ = nc_vector_shape_func(x+h, y, b, cntrd)
    ψ₁₋ = nc_vector_shape_func(x-h, y, b, cntrd)
    ψ₂₊ = nc_vector_shape_func(x, y+h, b, cntrd)
    ψ₂₋ = nc_vector_shape_func(x, y-h, b, cntrd)
    for j = 1:6
	∇ψ[:,1,j] = ( ψ₁₊[:,j] - ψ₁₋[:,j] ) / (2h)
	∇ψ[:,2,j] = ( ψ₂₊[:,j] - ψ₂₋[:,j] ) / (2h)
    end
end

function nc_nodes(a::Matrix{T}) where T <: AbstractFloat
    coord = Matrix{T}(undef, 2, 6)
    for k = 1:3
	k₊ = SUCC[k]
	coord[:,k] .= a[:,k]
	coord[:,k+3] .= ( a[:,k] + a[:,k₊] ) / 2
    end
    return coord
end

function contraction(A::Matrix{T}, B::Matrix{T}) where T <: Number
    return dot(A[:], B[:])
end

function conforming(λ, μ, coord, level, h)
    b, area = barycentric(coord)
    cntrd = centroid(coord)
    x, y = quadrature_points(level, coord)
    εψ = zeros(2, 2, 6)
    A = zeros(6, 6)
    B = zeros(6, 6)
    for k in eachindex(x)
	shape_strain!(εψ, x[k], y[k], b, cntrd, h)
        for j = 1:6
	    div_ψⱼ = εψ[1,1,j] + εψ[2,2,j]
	    for i = 1:6
	        div_ψᵢ = εψ[1,1,i] + εψ[2,2,i]
		A[i,j] += λ(x[k], y[k]) * div_ψⱼ * div_ψᵢ
		B[i,j] += 2μ(x[k], y[k]) * contraction(εψ[:,:,j], εψ[:,:,i])
	    end
	end
    end
    A *= area / length(x)
    B *= area / length(x)
    return A, B
end

function non_conforming(μ::Function, μ_plus_λ::Function, ∇μ::Function, 
	coord::Matrix{Float64}, level::Int64, h::Float64)
    b, area = barycentric(coord)
    cntrd = centroid(coord)
    x, y = quadrature_points(level, coord)
    A = zeros(6, 6)
    B = zeros(6, 6)
    D = zeros(6, 6)
    ∇ψ = zeros(2, 2, 6)
    for k in eachindex(x)
	nc_shape_gradient!(∇ψ, x[k], y[k], b, cntrd, h)
	ψ = nc_vector_shape_func(x[k], y[k], b, cntrd)
	v = ∇μ(x[k], y[k])
	for j = 1:6
	    div_ψⱼ = ∇ψ[1,1,j] + ∇ψ[2,2,j]
	    for i = 1:6
	        div_ψᵢ = ∇ψ[1,1,i] + ∇ψ[2,2,i]
		A[i,j] += μ_plus_λ(x[k], y[k]) * div_ψⱼ * div_ψᵢ
		B[i,j] += μ(x[k], y[k]) * contraction(∇ψ[:,:,j], ∇ψ[:,:,i])
		D[i,j] += (  v[2] * 
			     ( ∇ψ[1,1,j] * ψ[2,i] + ψ[2,j] * ∇ψ[1,1,i] )
			   - v[1] * 
			     ( ∇ψ[1,2,j] * ψ[2,i] + ψ[2,j] * ∇ψ[1,2,i] ) )
	    end
	end
    end
    A *= area / length(x)
    B *= area / length(x)
    D *= area / length(x)
    return A, B, D
end

function load_vector!(fvec::Vector{Float64}, f::Function, 
	coord::Matrix{Float64}, level::Int64)
    b, area = barycentric(coord)
    cntrd = centroid(coord)
    x, y = quadrature_points(level, coord)
    fill!(fvec, 0.0)
    for k in eachindex(x)
	ψ = vector_shape_func(x[k], y[k], b, cntrd)
	fk = f(x[k], y[k])
	for j = 1:6
	    fvec[j] += dot(fk, ψ[:,j])
	end
    end
    fvec .*= area / length(x)
end

function nc_load_vector!(fvec::Vector{Float64}, f::Function, 
	coord::Matrix{Float64}, level::Int64)
    b, area = barycentric(coord)
    cntrd = centroid(coord)
    x, y = quadrature_points(level, coord)
    fill!(fvec, 0.0)
    for k in eachindex(x)
	ψ = nc_vector_shape_func(x[k], y[k], b, cntrd)
	fk = f(x[k], y[k])
	for j = 1:6
	    fvec[j] += dot(fk, ψ[:,j])
	end
    end
    fvec .*= area / length(x)
end

let
    λ(x, y) = 2.0 + x + y/2
    μ(x, y) = 1.0 + x/3 + y
    μ_plus_λ(x, y) = μ(x, y) + λ(x, y)
    ∇μ(x, y) = SA[1/3, 1]
    f(x, y) = SA[ 1.0 - x + 3y, 2.0 + 2x - y/2 ]

    coord = [   π  -2.0  -1.0
              1.0   2.0  -2.0 ]

    A1 = Matrix{Float64}(undef, 6, 6)
    B1 = similar(A1)

    ∫∫λ_div_u_div_v!(A1, coord, λ)
    ∫∫2μ_εu_εv!(B1, coord, μ)
    b, area = barycentric(coord)
    cntrd = centroid(coord)

    h = 1e-3
    level = 1
    x, y = 1.5, 1.0
    εψ1 = Array{Float64}(undef, 2, 2, 6)
    εψ2 = similar(εψ1)
    shape_strain!(εψ1, x, y, b, cntrd, h)
    for j = 1:3
	εψ2[1,1,j] = b[1,j]
	εψ2[1,2,j] = εψ2[2,1,j] = b[2,j] / 2
	εψ2[2,2,j] = 0
	εψ2[1,1,j+3] = 0
	εψ2[1,2,j+3] = εψ2[2,1,j+3] = b[1,j] / 2
	εψ2[2,2,j+3] = b[2,j]
    end

    @test εψ1 ≈ εψ2

    A2, B2 = conforming(λ, μ, coord, level, h)

    @test A1 ≈ A2
    @test issymmetric(A1)

    @test B1 ≈ B2
    @test issymmetric(B1)

    f1 = Vector{Float64}(undef, 6)
    f2 = similar(f1)
    ∫∫f_dot_v!(f1, coord, f)
    load_vector!(f2, f, coord, level)

    @test f1 ≈ f2

    nc_coord = nc_nodes(coord)
    A3 = Matrix{Float64}(undef, 6, 6)
    B3 = similar(A3)
    D3 = similar(A3)
    NCE.∫∫a_div_u_div_v!(A3, nc_coord, μ_plus_λ)
    NCE.∫∫μ_∇u_colon_∇v!(B3, nc_coord, μ)
    NCE.correction!(D3, nc_coord, ∇μ)

    ∇ψ1 = Array{Float64}(undef, 2, 2, 6)
    ∇ψ2 = similar(∇ψ1)
    nc_shape_gradient!(∇ψ1, x, y, b, cntrd, h)
    for j = 1:3
	j₋ = PRED[j]
	∇ψ2[1,1,j] = -2b[1,j₋]
	∇ψ2[1,2,j] = -2b[2,j₋]
	∇ψ2[2,:,j] .= 0
	∇ψ2[1,:,j+3] .= 0
	∇ψ2[2,1,j+3] = ∇ψ2[1,1,j]
	∇ψ2[2,2,j+3] = ∇ψ2[1,2,j]
    end

    @test ∇ψ1 ≈ ∇ψ2

    A4, B4, D4 = non_conforming(μ, μ_plus_λ, ∇μ, coord, level, h)

    @test A3 ≈ A4
    @test issymmetric(A3)

    @test B3 ≈ B4
    @test issymmetric(B3)

    @test D3 ≈ D4
    @test issymmetric(D3)

    f3 = Vector{Float64}(undef, 6)
    f4 = similar(f1)
    NCE.∫∫f_dot_v!(f3, nc_coord, f)
    nc_load_vector!(f4, f, coord, level)

    @test f3 ≈ f4
end
