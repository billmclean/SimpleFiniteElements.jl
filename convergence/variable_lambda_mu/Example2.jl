module Example2

import StaticArrays: SA

Λ = 2.0
ϵ = 1.0
scale = 10.0  # Increase (descrease) to make arrows smaller (larger).

function λ(x, y)
    return Λ * ( 2 + sin(2x) )
end

function μ(x, y)
    return 1 + ϵ * ( x + y )
end

function f(x, y)
    w1 = ( ( 2 * sin(x) * sin(y) - cos(x+y) ) * μ(x, y) / Λ
		 - cos(x+y) * ( sin(2x) + 2 ) )
    w2 = cos(2x) - cos(2y)
    w3 = 2 * sin(2x) * sin(2y)
    f1 = ( sin(x-y) - sin(3x+y) + w1 + 4 * sin(2y) * ( 2*cos(2x) - 1 ) * μ(x, y)
	  - ϵ * ( 2 * cos(x) * sin(y) + sin(x+y) ) / Λ + 2ϵ * ( w3 - w2 ) )
    f2 = ( w1 - 4 * sin(2x) * ( 2 * cos(2y) - 1 ) * μ(x, y)
	  - ϵ * ( 2 * sin(x) * cos(y) + sin(x+y) ) / Λ
	  - 2ϵ * ( w2 + w3 ) )
    return SA[f1, f2]
end

function u(x, y)
    v = sin(x) * sin(y) / Λ
    u1 = ( cos(2x) - 1 ) * sin(2y) + v
    u2 = ( 1 - cos(2y) ) * sin(2x) + v
    return SA[u1, u2]
end

msg = """Example 2. 
λ(x, y) = Λ ( sin(2x) + 2 ) for Λ = $Λ.
μ(x, y) = 1 + ϵ ( x + y ) for ϵ = $ϵ."""

end # module
