Point(1) = { 0,   0, 0 };
Point(2) = { 48, 44, 0 };
Point(3) = { 48, 60, 0 };
Point(4) = { 0,  44, 0 };
Line(1) = { 1, 2 };
Line(2) = { 2, 3 };
Line(3) = { 3, 4 };
Line(4) = { 4, 1 };
Line Loop(1) = { 1, 2, 3, 4 };
Plane Surface(1) = { 1 };

Physical Surface("Omega") = { 1 };

Physical Line("Bottom") = { 1 };
Physical Line("Right")  = { 2 };
Physical Line("Top")    = { 3 };
Physical Line("Left")   = { 4 };
