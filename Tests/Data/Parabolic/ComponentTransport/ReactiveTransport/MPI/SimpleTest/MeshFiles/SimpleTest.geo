//SetFactory("OpenCASCADE");
Point(1) = {0, 0, 0, 1};
Point(2) = {1, 0, 0, 1};
Point(3) = {1, 1, 0, 1};
Point(4) = {0, 1, 0, 1};

Line(1) = {1, 2}; Transfinite Line {1} = 11 Using Progression 1;
Line(2) = {2, 3}; Transfinite Line {2} = 11 Using Progression 1;
Line(3) = {3, 4}; Transfinite Line {3} = 11 Using Progression 1;
Line(4) = {4, 1}; Transfinite Line {4} = 11 Using Progression 1;

// Surface
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Structured mesh
Transfinite Surface {1};
Recombine Surface {1};

// Physical group
Physical Surface(1) = {1};
