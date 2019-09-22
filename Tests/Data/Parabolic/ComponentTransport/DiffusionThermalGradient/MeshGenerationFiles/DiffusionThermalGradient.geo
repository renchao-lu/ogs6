// Grid number
GridNum = 41;

// Points
Point(1) = {0, 0, 0, 1.0};
Point(2) = {0.2, 0, 0, 1.0};
Point(3) = {0.2, 10, 0, 1.0};
Point(4) = {0, 10, 0, 1.0};
Point(5) = {10, 0, 0, 1.0};
Point(6) = {10, 10, 0, 1.0};
Point(7) = {10, 15, 0, 1.0};
Point(8) = {0.2, 15, 0, 1.0};

// Lines
Line(1) = {1, 2}; Transfinite Line {1} = 5 Using Progression 1;
Line(2) = {2, 3}; Transfinite Line {2} = 201 Using Progression 1;
Line(3) = {3, 4}; Transfinite Line {3} = 5 Using Progression 1;
Line(4) = {4, 1}; Transfinite Line {4} = 201 Using Progression 1;
Line(5) = {2, 5}; Transfinite Line {5} = 197 Using Progression 1;
Line(6) = {5, 6}; Transfinite Line {6} = 201 Using Progression 1;
Line(7) = {6, 3}; Transfinite Line {7} = 197 Using Progression 1;
Line(8) = {6, 7}; Transfinite Line {8} = 101 Using Progression 1;
Line(9) = {7, 8}; Transfinite Line {9} = 197 Using Progression 1;
Line(10) = {8, 3}; Transfinite Line {10} = 101 Using Progression 1;

// Surfaces
Line Loop(11) = {1, 2, 3, 4};
Plane Surface(12) = {11};
Line Loop(13) = {5, 6, 7, -2};
Plane Surface(14) = {13};
Line Loop(15) = {8, 9, 10, -7};
Plane Surface(16) = {15};

// For structured mapping
Transfinite Surface {12};
Transfinite Surface {14};
Transfinite Surface {16};
Recombine Surface {12};
Recombine Surface {14};
Recombine Surface {16};

// Physical groups
Physical Surface(47) = {12, 14, 16};
