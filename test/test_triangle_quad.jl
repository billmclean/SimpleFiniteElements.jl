import SimpleFiniteElements.Utils: barycentric, triangle_quad
import LinearAlgebra: dot

let
    a = [ 2.0  -2.0  -1.0
          1.0   2.0  -2.0 ]

    b, area = barycentric(a)
    f(x, y) = x * y

    Q = triangle_quad(a, f)
    I = area * Q

    diag = a[1,1] * a[2,1] + a[1,2] * a[2,2] + a[1,3] * a[2,3] 
    off_diag = (                    a[1,1] * a[2,2] + a[1,1] * a[2,3]
    	        + a[1,2] * a[2,1]                   + a[1,2] * a[2,3]
	        + a[1,3] * a[2,1] + a[1,3] * a[2,2]                   )

    exact = (area/6) * diag + (area/12) * off_diag

    @test abs(I - exact) < eps(I)
end
