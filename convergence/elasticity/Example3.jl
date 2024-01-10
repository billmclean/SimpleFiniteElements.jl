module Example3

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
a = 1.5
b = 1.8

msg = """
Example 3: elasticity equation with constant coefficients, non-zero load.
Spatial domain defined in $path.
Lame parameters λ = $λ, μ = $μ.
Exact solution
    u = [ x^2 + xy - a sin(y), xy - y^2 + b exp(x) ]
with a = $a and b = $b."""

exact_u(x, y) = [ x^2 + x*y - a * sin(y) 
		  x*y - y^2 + b * exp(x) ]
g = exact_u 

f(x, y) = [ -(3λ + 5μ + a * μ * sin(y)), λ + 3μ - b * μ * exp(x) ]
essential_bcs = [("Gamma", g)]

end # module
