cl = .05;
L = 1.0;
alpha = 0.4;

Point(1) = {0, 0, 0, cl};
Point(2) = {L, 0, 0, cl};
Point(3) = {L, L, 0, cl};
Point(4) = {0, L, 0, cl};
Point(5) = {0, L, alpha*L, cl};
Point(6) = {0, 0, alpha*L, cl};
Point(10) = {L, 0, alpha*L, cl};
Point(14) = {L, L, alpha*L, cl};

Line(1) = {4, 1};
Line(2) = {1, 2};
Line(3) = {2, 3};
Line(4) = {3, 4};
Line(8) = {5, 6};
Line(9) = {6, 10};
Line(10) = {10, 14};
Line(11) = {14, 5};
Line(12) = {6, 1};
Line(13) = {4, 5};
Line(14) = {14, 3};
Line(15) = {10, 2};

Line Loop(6) = {1, 2, 3, 4};
Plane Surface(6) = {6};

Line Loop(7) = {8, 9, 10, 11};
Plane Surface(7) = {7};

Line Loop(17) = {8, 12, -1, 13};
Plane Surface(17) = {17};

Line Loop(19) = {13, -11, 14, 4};
Plane Surface(19) = {19};

Line Loop(21) = {14, -3, -15, 10};
Plane Surface(21) = {21};

Line Loop(23) = {15, -2, -12, 9};
Plane Surface(23) = {23};

Surface Loop(1)={6,7,17,19,21,23};
Volume(1)={1};

Physical Surface(10) = {23};					// Poco Injetor
Physical Surface(51) = {19};					// Poco Produtor
Physical Surface(2000) = {6,7,17,21};			// Faces de Contorno
Physical Volume(3000) = {1};					// Define meio poroso
