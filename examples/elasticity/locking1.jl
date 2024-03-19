using SimpleFiniteElements
import SimpleFiniteElements.Elasticity: ∫∫f_dot_v!, ∫∫λ_div_u_div_v!, 
                                        ∫∫2μ_εu_εv!, elasticity_soln, error_norms
import Printf: @printf
import StaticArrays: SA

function exact_u(x, y, λ)
    common_term = sin(x) * sin(y) / λ
    u1 = ( cos(2x) - 1 ) * sin(2y) + common_term
    u2 = ( 1 - cos(2y) ) * sin(2x) + common_term
    return SA[u1, u2]
end

function ∇u(x, y, λ)
    sx, cx = sincos(x)
    sy, cy = sincos(y)
    s2x, c2x = sincos(2x)
    s2y, c2y = sincos(2y)
    s2x_s2y = s2x * s2y
    cx_sy = cx * sy
    sx_cy = sx * cy
    ∂₁u₁ = -2s2x_s2y        + cx_sy / λ
    ∂₂u₁ = 2(c2x - 1) * c2y + sx_cy / λ
    ∂₁u₂ = 2(1 - c2y) * c2x + cx_sy / λ
    ∂₂u₂ =  2s2x_s2y        + sx_cy / λ
    return SA[ ∂₁u₁  ∂₁u₂
               ∂₂u₁  ∂₂u₂ ]
end

function ϕ(x, y, μ)
    sx, cx = sincos(2x)
    sy, cy = sincos(2y)
    return 4μ * sy * (2cx - 1)
end

function f(x, y, λ, μ)
    ψ1 = μ * ( 2sin(x) * sin(y) - cos(x+y) ) 
    ψ3 = λ * cos(x+y)
    f1 =  ϕ(x, y, μ) - ψ3 + ψ1 / λ
    f2 = -ϕ(y, x, μ) - ψ3 + ψ1 / λ
    return SA[f1, f2]
end

path = joinpath("..", "..", "spatial_domains", "square.geo")
gmodel = GeometryModel(path)

λ = 1.0
μ = 2.0
gD = SA[0.0, 0.0]
essential_bcs = [("Top", gD), ("Bottom", gD), ("Left", gD), ("Right", gD)]
bilinear_forms = Dict("Omega" => [(∫∫λ_div_u_div_v!, λ),
                                  (∫∫2μ_εu_εv!, μ)])
linear_funcs = Dict("Omega" => (∫∫f_dot_v!, f, λ, μ))

hmax = 0.2
nrows = 4
mesh = FEMesh(gmodel, hmax, order=1, save_msh_file=false,
              refinements=nrows-1, verbosity=2)

L2err = Vector{Float64}(undef, nrows)
H1err = similar(L2err)

quadrature_level = 2
@printf("Conforming elements, λ = %g, μ = %g, quadrature_level = %g.\n", λ, μ, quadrature_level)
@printf("%6s& %6s& %8s& %5s& %8s& %5s\\\\\n\n",
        "h", "DoF", "L2", "rate", "H1", "rate")
for row in eachindex(L2err)
    start = time()
    dof = DegreesOfFreedom(mesh[row], essential_bcs)
    u1h, u2h = elasticity_soln(dof, bilinear_forms, linear_funcs)
    L2err[row], H1err[row] = error_norms(u1h, u2h, exact_u, ∇u,
                                         dof, quadrature_level, λ)
    h = max_elt_diameter(mesh[row])
    num_dof = 2 * dof.num_free
    elapsed = time() - start
    if row == 1
        @printf("%6.3f& %6d& %8.2e& %5s& %8.2e& %5s\\\\\n",
                h, num_dof, L2err[row], "",    H1err[row], "")
    else
        L2rate = log2(L2err[row-1] / L2err[row])
        H1rate = log2(H1err[row-1] / H1err[row])
        @printf("%6.3f& %6d& %8.2e& %5.3f& %8.2e& %5.3f\\\\\n",
                h, num_dof, L2err[row], L2rate,    H1err[row], H1rate)
    end
end

