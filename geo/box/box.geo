// gmsh -2 -smooth 1000 -clmin 0.75 -clmax 2 disk.geo -o disk.msh

// control points
Point(1) = {-1, -1, 0, 0.5};
Point(2) = {+1, -1, 0, 0.5};
Point(3) = {+1, +1, 0, 0.5};
Point(4) = {-1, +1, 0, 0.5};

// North
Line(5) = {1, 2};
Line(6) = {2, 3};
Line(7) = {3, 4};
Line(8) = {4, 1};

Line Loop(9) = {5, 6, 7, 8};
Physical Line("boundary", 1) = {5, 6, 7, 8};

Plane Surface(10) = {9};
Mesh.CharacteristicLengthMin = 10;
Mesh.CharacteristicLengthMax = 10;

// 2D mesh algorithm (1=MeshAdapt, 2=Automatic, 5=Delaunay, 6=Frontal, 7=bamg, 8=delquad)
Mesh.Algorithm = 1;

Physical Surface("elm", 999) = {10};

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
RefineMesh;
Save StrCat(StrPrefix(General.FileName), "_r6.msh");
Exit;
