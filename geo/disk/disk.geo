// Circle radius
R = 1.0;

// Inner box size
L = Sqrt(2)/2;

// control points
Point(1) = { 0,  0, 0, 0.5};
Point(2) = {+L, +L, 0, 0.5};
Point(3) = {-L, +L, 0, 0.5};
Point(4) = {-L, -L, 0, 0.5};
Point(5) = {+L, -L, 0, 0.5};

// North
Circle(1) = {2, 1, 3};
Line(2) = {3, 2};
Line Loop(3) = {1, 2};

// West
Circle(4) = {3, 1, 4};
Line(5) = {4, 3};
Line Loop(6) = {4, 5};

// South
Circle(7) = {4, 1, 5};
Line(8) = {5, 4};
Line Loop(9) = {7, 8};

// East
Circle(10) = {5, 1, 2};
Line(11) = {2, 5};
Line Loop(12) = {11, 10};

Line Loop(13) = {2,5,8,11};

Plane Surface(14) = {3};
Plane Surface(15) = {6};
Plane Surface(16) = {9};
Plane Surface(17) = {12};
Plane Surface(18) = {13};

Physical Line("boundary", 1) = {1, 4, 7, 10};

// 2D mesh algorithm (1=MeshAdapt, 2=Automatic, 5=Delaunay, 6=Frontal, 7=bamg, 8=delquad)
Mesh.Algorithm = 1;
Mesh.CharacteristicLengthMin = 0.5;
Mesh.CharacteristicLengthMax = 0.5;

// Mesh recombination algorithm (0=standard, 1=blossom)
Mesh.RecombinationAlgorithm = 1;

Physical Surface("elm", 999) = {14, 15, 16, 17, 18};

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
