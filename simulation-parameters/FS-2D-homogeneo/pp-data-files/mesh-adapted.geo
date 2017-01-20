cl = 0.1;

Point(1) = {1, 0, 0, cl};
Point(2) = {1, 1, 0, cl};
Point(3) = {0, 1, 0, cl};
Point(4) = {0, 0, 0, cl};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1) = {2, 3, 4, 1};
Plane Surface(1) = {1};

Physical Line(2000) = {2,4};
Physical Line(51) = {1};
Physical Line(10) = {3};

Physical Surface(3300) = {1};
