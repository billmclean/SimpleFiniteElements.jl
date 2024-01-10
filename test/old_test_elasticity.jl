import SimpleFiniteElements.Elasticity: ∫∫λ_div_u_div_v!, ∫∫2μ_εu_εv!, 
                                        ∫∫f_dot_v!
import LinearAlgebra: norm
import StaticArrays: SA

coord = [ 0.0  2.0  1.0
          0.0  1.0  3.0 ]

A1 = zeros(6, 6)
λ = 1.0
∫∫λ_div_u_div_v!(A1, coord, λ)

A2 = zeros(6, 6)
∫∫λ_div_u_div_v!(A2, coord, (x,y) -> λ)

@test maximum(abs, A1-A2) < eps(norm(A1, Inf))

B1 = zeros(6, 6)
μ = 1.0
∫∫2μ_εu_εv!(B1, coord, μ)
B2 = zeros(6, 6)
∫∫2μ_εu_εv!(B2, coord, (x, y) -> μ)

@test maximum(abs, B1-B2) < eps(norm(B1, Inf))

f1 = SA[1.5, -0.75]
v1 = zeros(6)
∫∫f_dot_v!(v1, coord, f1)

v2 = zeros(6)
∫∫f_dot_v!(v2, coord, (x, y) -> f1)
