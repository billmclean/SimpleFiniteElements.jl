module Example1

using SimpleFiniteElements
import StaticArrays: SA
import SimpleFiniteElements.Elasticity: fundamental_soln
import Printf: @printf

export label, gmodel, λ, μ, exact_u, f, g, essential_bcs

path = joinpath("..", "..", "spatial_domains", "house.geo")
gmodel = GeometryModel(path)
λ = 2.0
μ = 1.0
μ_plus_λ = μ + λ
∇μ = SA[0.0, 0.0]

source_point = SA[0.4, 1.2]
source_vector = SA[2.0, -4.0]

msg = """
Example 1: elasticity equation with constant coefficients, zero load.
Spatial domain defined in $path.
Lame parameters λ = $λ, μ = $μ.
Source point = $source_point.
Source vector = $source_vector."""

g(x, y) = fundamental_soln(SA[x, y], source_point, λ, μ) * source_vector
exact_u = g

f(x, y) = SA[0.0, 0.0]
essential_bcs = [("Gamma", g)]

end # module
