import SimpleFiniteElements.NonConformingPoisson: ∫∫a_∇u_dot_∇v!

h = 2.0

coord = zeros(2, 6)
coord[:,1] = [0.0, 0.0]
coord[:,2] = [h, 0.0]
coord[:,3] = [0.0, h]
coord[:,4] = [h/2, 0.0]
coord[:,5] = [h/2, h/2]
coord[:,6] = [0.0, h/2]

A = Matrix{Float64}(undef, 3, 3)
coef = 2.0
∫∫a_∇u_dot_∇v!(A, coord, coef)

println("A = ")
display(A)

