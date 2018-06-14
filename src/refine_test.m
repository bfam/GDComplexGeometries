addpath ../geo/inclusion
addpath ../nodal-dg/Codes1.1/Codes1D
addpath ../nodal-dg/Codes1.1/Codes2D
addpath ../nodal-dg/Codes1.1/ServiceRoutines

clear
Globals2D;

N = 4;
refinelevel = 1;
mesh = 'inclusion_r0.msh';

[Nv, VX, VY, K, EToV, BCType] = MeshReaderGmshBC2D(mesh);
StartUp2D;

xout = RefineUniformField2D(refinelevel, x);
yout = RefineUniformField2D(refinelevel, y);

save('test.mat', 'N', 'mesh', 'refinelevel', 'xout', 'yout');

clear
Globals2D;
load('test.mat');

[Nv, VX, VY, K, EToV, BCType] = MeshReaderGmshBC2D(mesh);
empty = @(x) 0;
RefineUniform2D(refinelevel, empty);

StartUp2D;

MassMatrix = invV'*invV;

errx = x - xout;
erry = y - yout;

MMerrx = MassMatrix*(J.*errx);
MMerry = MassMatrix*(J.*erry);

L2errx = sqrt(errx(:)'*MMerrx(:))
L2erry = sqrt(erry(:)'*MMerry(:))
