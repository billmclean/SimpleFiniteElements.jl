module Example2

using SimpleFiniteElements
import StaticArrays: SA
import Printf: @printf

export gmodel, Λ, λ, μ, exact_u, f, g, essential_bcs, msg

path = joinpath("..", "..", "spatial_domains", "square.geo")
gmodel = GeometryModel(path)

Λ = 1.0
α = 0.5
τ(x, y) = 1 + α * sin(2x)
λ(x, y) = Λ * τ(x, y)

β = 1.0
μ(x, y) = 1 + β * (x + y)
∇μ = SA[β, β]

μ_plus_λ(x, y) = μ(x, y) + λ(x, y)

msg = """
Example 2: elasticity equation with variable coefficients.
Spatial domain defined in $path.
λ(x, y) = Λ ( 1 + α sin(2πx) ) with Λ = $Λ and α = $α
μ(x, y) = 1 + β(x + y)         with β = $β"""

function exact_u(x, y)
    common_term = sin(x) * sin(y) / Λ
    u1 = ( cos(2x) - 1 ) * sin(2y) + common_term
    u2 = ( 1 - cos(2y) ) * sin(2x) + common_term
    return SA[u1, u2]
end

function ϕ(x, y)
    sx = sin(2x)
    sy = sin(2y)
    cx = cos(2x)
    cy = cos(2y)
    return 2β * ( 2sx * sy - cx + cy ) + 4μ(x, y) * sy * (2cx - 1)
end

function f(x, y)
    ψ1 = μ(x, y) * ( 2sin(x) * sin(y) - cos(x+y) ) - β * sin(x+y)
    ψ2(x, y) = 2β * cos(x) * sin(y)
    ψ3 = τ(x, y) * cos(x+y)
    f1 =  ϕ(x, y) - ψ3 - 2α * cos(2x) * sin(x+y) + (ψ1 -  ψ2(x, y)) / Λ
    f2 = -ϕ(y, x) - ψ3 + (ψ1 - ψ2(y, x)) / Λ
    return SA[f1, f2]
end

g(x, y) = SA[0.0, 0.0]
#g = exact_u

essential_bcs = [("Top", g), ("Bottom", g), ("Left", g), ("Right", g)]

end # module
