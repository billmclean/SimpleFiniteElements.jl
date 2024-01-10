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

Line(1) = { 8, 1 };
Line(2) = { 1, 2 };
Circle(3) = { 2, 9, 3 };
Line(4) = { 3, 4 };
Line(5) = { 4, 5 };
Line(6) = { 5, 6 };
Circle(7) = { 6, 10, 7 };
Line(8) = { 7, 8 };

Line Loop(1) = { 1, 2, 3, 4, 5, 6, 7, 8 };
Plane Surface(1) = { 1 };

Physical Surface("Omega") = { 1 };
Physical Line("Fixed") = { 1 };
Physical Line("Traction") = { 5 };
Physical Line("NoTraction") = { 2, 3, 4, 6, 7, 8 };

