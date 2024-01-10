Point(1) = {-2, -1, 0};
Point(2) = {-1, -1, 0};
Point(3) = { 1, -1, 0};
Point(4) = { 2, -1, 0};
Point(5) = { 2,  1, 0};
Point(6) = { 1,  1, 0};
Point(7) = {-1,  1, 0};
Point(8) = {-2,  1, 0};
R = 2;
Point(9)  = { 0, -0.5-R, 0 };
Point(10) = { 0,  0.5+R, 0 };

Line(1) = { 3, 4, 5, 6 };
Circle(2) = { 6, 10, 7 };
Line(3) = { 7, 8, 1, 2 };
Circle(4) = { 2, 9, 3 };

Line Loop(1) = { 1, 2, 3, 4 };
Plane Surface(1) = { 1 };

Physical Surface("Omega") = { 1 };
Physical Line("Gamma") = { 1, 2, 3, 4 };

