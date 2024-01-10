a = 2;
h = 1;
Point(1) = { -a, -a, 0 };
Point(2) = {  a, -a, 0 };
Point(3) = {  a,  0, 0 };
Point(4) = {  0,  h, 0 };
Point(5) = { -a,  0, 0 };

Line(1) = { 1, 2 };
Line(2) = { 2, 3 };
Line(3) = { 3, 4 };
Line(4) = { 4, 5 };
Line(5) = { 5, 1 };

Line Loop(1) = { 1, 2, 3, 4, 5 };
Plane Surface(1) = { 1 };

Physical Surface("Omega") = { 1 };
Physical Line("Gamma") = { 1, 2, 3, 4, 5 };
