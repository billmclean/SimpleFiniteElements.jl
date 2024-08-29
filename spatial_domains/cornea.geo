R = 3;
r = 2;
alpha = Pi/6;

Point(1) = { 0, 0, 0 };
Point(2) = {-r * Sin(alpha), r * Cos(alpha), 0 };
Point(3) = { r * Sin(alpha), r * Cos(alpha), 0 };
Point(4) = { R * Sin(alpha), R * Cos(alpha), 0 };
Point(5) = {-R * Sin(alpha), R * Cos(alpha), 0 };

Circle(1) = { 2, 1, 3 };
Line(2) = { 3, 4 };
Circle(3) = { 4, 1, 5 };
Line(4) = { 5, 2 };

Line Loop(1) = { 1, 2, 3, 4 };
Plane Surface(1) = { 1 };

Physical Surface("Omega") = { 1 };
Physical Line("Gamma") = { 1, 2, 3, 4 };
