Point(1) = { -2.0,  -2.0,  0.0 };
Point(2) = {  2.0,  -2.0,  0.0 };
Point(3) = {  2.0,   2.0,  0.0 };
Point(4) = { -2.0,   2.0,  0.0 };
Point(5) = {  0.0,   0.0,  0.0 };
Point(6) = {  0.0,  -1.0,  0.0 };
Point(7) = {  1.0,   0.0,  0.0 };
Point(8) = {  0.0,   1.0,  0.0 };
Point(9) = { -1.0,   0.0,  0.0 };

Line(1) = { 1, 2 };
Line(2) = { 2, 3 };
Line(3) = { 3, 4 };
Line(4) = { 4, 1 };
Circle(5) = {6, 5, 7 };
Circle(6) = {7, 5, 8 };
Circle(7) = {8, 5, 9 };
Circle(8) = {9, 5, 6 };

Line Loop(1) = { 1, 2, 3, 4 };
Line Loop(2) = { 5, 6, 7, 8 };
Plane Surface(1) = { 1, -2 };
Plane Surface(2) = { 2 };

Physical Line("Outer Bdry") = { 1, 2, 3, 4 };
Physical Line("Inner Bdry") = { -8, -7, -6, -5 };
Physical Surface("Outer Domain") = { 1 };
Physical Surface("Inner Domain") = { 2 };
