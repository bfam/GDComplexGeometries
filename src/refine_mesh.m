addpath ../geo/inclusion
addpath ../nodal-dg/Codes1.1/Codes1D
addpath ../nodal-dg/Codes1.1/Codes2D
addpath ../nodal-dg/Codes1.1/ServiceRoutines

clear

Globals2D;

mesh = 'inclusion_r0.msh';

[Nv, VX, VY, K, EToV, BCType] = MeshReaderGmshBC2D(mesh);

RefineUniform2D(1, @MakeInclusionsVerts2D);

N = 1;
StartUp2D;

WriteVTK2D('mesh.vtk',1,{});
