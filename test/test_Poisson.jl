import SimpleFiniteElements.Poisson: ∫∫a_∇u_dot_∇v!, ∫∫c_u_v!, 
                                     ∫c_u_v!, ∫∫f_v!, ∫g_v!

let
    coord = [ 0.0  2.0  1.0
              0.0  1.0  3.0 ]

    A1 = zeros(3, 3)
    ∫∫a_∇u_dot_∇v!(A1, coord, 2.0)
    A2 = zeros(3, 3)
    ∫∫a_∇u_dot_∇v!(A2, coord, (x,y)->2.0)
    @test maximum(abs, A1-A2) < eps(Float64)

    M1 = zeros(3, 3)
    ∫∫c_u_v!(M1, coord, 2.0)
    M2 = zeros(3, 3)
    ∫∫c_u_v!(M2, coord, (x,y)->2.0)
    @test maximum(abs, M1-M2) < eps(Float64)

    coordE = [ 0.0  2.0
               0.0  1.0 ]
    ME1 = zeros(2, 2)
    ∫c_u_v!(ME1, coordE, 3.0)
    ME2 = zeros(2, 2)
    ∫c_u_v!(ME2, coordE, (x,y)->3.0)
    @test maximum(abs, ME1-ME2) < eps(Float64)

    v1 = zeros(3)
    ∫∫f_v!(v1, coord, 2.0)
    v2 = zeros(3)
    ∫∫f_v!(v2, coord, (x,y)->2.0)
    @test maximum(abs, v1-v2) < eps(Float64)

    vE1 = zeros(2)
    ∫g_v!(vE1, coordE, 2.0)
    vE2 = zeros(2)
    ∫g_v!(vE2, coordE, (x,y)->2.0)
    @test maximum(abs, vE1-vE2) < eps(Float64)
end
