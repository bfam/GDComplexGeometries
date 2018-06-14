
// outside box
r1 = 5;

Point(1) = {-r1,-r1,0};
Point(2) = { r1,-r1,0};
Point(3) = { r1, r1,0};
Point(4) = {-r1, r1,0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Physical Line("south", 3) = {1};
Physical Line("west", 1) = {4};
Physical Line("north", 4) = {3};
Physical Line("east", 2) = {2};

// cylinder 2
r2 = 2;
x2 = 2;
y2 = 2;

Point(5) = {x2+ 0,y2+ 0,0};

Point(6) = {x2+ 0,y2-r2,0};
Point(7) = {x2+ 0,y2+r2,0};
Point(8) = {x2+r2,y2+ 0,0};
Point(9) = {x2-r2,y2+ 0,0};

Circle(5) = {8, 5, 7};
Circle(6) = {7, 5, 9};
Circle(7) = {9, 5, 6};
Circle(8) = {6, 5, 8};

Physical Line("cyl2", 5) = {5,6,7,8};

// cylinder 3
r3 =  2.5;
x3 = -2;
y3 = -2;

Point(10) = {x3+ 0,y3+ 0,0};

Point(11) = {x3+ 0,y3-r3,0};
Point(12) = {x3+ 0,y3+r3,0};
Point(13) = {x3+r3,y3+ 0,0};
Point(14) = {x3-r3,y3+ 0,0};

Circle( 9) = {13, 10, 12};
Circle(10) = {12, 10, 14};
Circle(11) = {14, 10, 11};
Circle(12) = {11, 10, 13};

Physical Line("cyl3", 6) = {9,10,11,12};

// cylinder 4
r4 = 1.5;
x4 = -3;
y4 = 2.5;

Point(15) = {x4+ 0,y4+ 0,0};

Point(16) = {x4+ 0,y4-r4,0};
Point(17) = {x4+ 0,y4+r4,0};
Point(18) = {x4+r4,y4+ 0,0};
Point(19) = {x4-r4,y4+ 0,0};

Circle(13) = {18, 15, 17};
Circle(14) = {17, 15, 19};
Circle(15) = {19, 15, 16};
Circle(16) = {16, 15, 18};

Physical Line("cyl4", 7) = {13,14,15,16};

// cylinder 5
r5 = 1;
x5 = 2.5;
y5 = -3;

Point(20) = {x5+ 0,y5+ 0,0};

Point(21) = {x5+ 0,y5-r5,0};
Point(22) = {x5+ 0,y5+r5,0};
Point(23) = {x5+r5,y5+ 0,0};
Point(24) = {x5-r5,y5+ 0,0};

Circle(17) = {23, 20, 22};
Circle(18) = {22, 20, 24};
Circle(19) = {24, 20, 21};
Circle(20) = {21, 20, 23};

Physical Line("cyl5", 8) = {17,18,19,20};

Line Loop(21) = { 1, 2, 3, 4}; // outside box 1
Line Loop(22) = { 5, 6, 7, 8}; // cylinder 2
Line Loop(23) = { 9,10,11,12}; // cylinder 3
Line Loop(24) = {13,14,15,16}; // cylinder 4
Line Loop(25) = {17,18,19,20}; // cylinder 5

Plane Surface(16) = {21, 22, 23, 24, 25};
// Recombine Surface {11};

// 2D mesh algorithm (1=MeshAdapt, 2=Automatic, 5=Delaunay, 6=Frontal, 7=bamg, 8=delquad)
Mesh.Algorithm = 7;

// Mesh recombination algorithm (0=standard, 1=blossom)
Mesh.RecombinationAlgorithm = 1;

// Mesh.RecombineAll = 1;

// Characteristic Length {17, 19, 16, 18} = 0.1;
//+
Physical Surface("elm", 999) = {16};

Mesh 2;

Save StrCat(StrPrefix(General.FileName), "_r0.msh");
RefineMesh;
Save StrCat(StrPrefix(General.FileName), "_r1.msh");
RefineMesh;
Save StrCat(StrPrefix(General.FileName), "_r2.msh");
RefineMesh;
Save StrCat(StrPrefix(General.FileName), "_r3.msh");
RefineMesh;
Save StrCat(StrPrefix(General.FileName), "_r4.msh");
RefineMesh;
Save StrCat(StrPrefix(General.FileName), "_r5.msh");
Exit;
