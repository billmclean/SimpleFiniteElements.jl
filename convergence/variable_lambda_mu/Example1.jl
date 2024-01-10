module Example1

import StaticArrays: SA

function λ(x, y, Λ)
    return Λ * ( 2 + sin(2x) )
end

function μ(x, y, α)
    return 1 + α * (x + y)
end

function ϕ(x, y, α)
    sx = sin(2x)
    sy = sin(2y)
    cx = cos(2x)
    cy = cos(2y)
    return 2α * ( 2sx * sy - cx + cy ) + 4μ(x, y, α) * sy * (2cx - 1)
end

function u(x, y)
    u1 = ( cos(2x) - 1 ) * sin(2y)
    u2 = ( 1 - cos(2y) ) * sin(2x)
    return SA[u1, u2]
end

function f(x, y, α, Λ)
    return SA[ϕ(x, y, α), -ϕ(y, x, α)]
end

msg = """Example 1. 
λ(x, y) = Λ ( 2 + sin(2πx) ) 
μ(x, y) = 1 + α(x + y)""" 

end # module
