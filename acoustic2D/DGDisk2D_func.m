% Curved, noncoforming simplicial DG Disk test
% Inputs:
%   p         :: polynomial order is 2*p+1
%   sr        :: level of refiment, e.g, N-> N0*2^sr
%   mesh_base :: mesh used is sprintf('%s_r%d.msh', mesh_base, sr)
%   root      :: (optional) which root for the solution to use (default: [])
%                NOTE: if isempty(root) then the codes does not solve the disk
%                      problem but instead computes the operator matrix
% Outputs:
%   T    :: time
%   err  :: energy
function [T,err] = DGDisk2D_func(p, sr, mesh_base, root)

fprintf('DGDisk2D_func(%3d, %3d, %s)\n', p, sr, mesh_base)

addpath ../src
addpath ../geo/disk
Globals2D_gddg
addpath([NODAL_DG_ROOT, '/Codes1.1/Codes1D'])
addpath([NODAL_DG_ROOT, '/Codes1.1/Codes2D'])
addpath([NODAL_DG_ROOT, '/Codes1.1/ServiceRoutines'])

Globals2D;

% order of the DG
N = 2*p+1;

% order of the mortar
Nm = 2*p+1;

% DG initialization
mesh = sprintf('%s_r%d.msh', mesh_base, sr);
[Nv, VX, VY, K, EToV, BCType] = MeshReaderGmshBC2D(mesh);
StartUp2D;
% BuildBCMaps2D
BuildBCMaps2D_gddg()
vmapDG = reshape([1:Np*K], size(x));

% setup the curved grid
[k, f] = find(BCType == 1);
if(~isempty(k))
  outfaces = [k,f];
  MakeCylinder2D(outfaces, 1, 0, 0);
end
intC = 2*(N+1) + 1;
OP.dg.cubature = CubatureVolumeMesh2D(intC);
OP.dg.cubature.WJI = OP.dg.cubature.W ./ (OP.dg.cubature.J).^2;

intG = 2*(N+1);
OP.dg.gauss = GaussFaceMesh2D(intG);
[OP.dg.gauss.mapB, OP.dg.gauss.mapMor] = ...
  BuildBCMaps2D_cubature(intG, BCType, OP.dg.gauss.mapM);

% Pretend the DG if the last GD block
OP.B{1}.isdg = true;
OP.B{1}.toB = [];
OP.B{1}.toF = [];
OP.B{1}.quad_order = Nfp-1;
OP.B{1}.vmap = vmapDG;

  % Create mortar grids for the DG block

  % Here we assume that the "mortar" is aligned with one of the coordinate
  % directions with the mapping that the DG faces match the GD faces.

  % Set up the global grid and initial conditions
  x1 = [x(:)];
  x2 = [y(:)];
  OP.x1 = x1; OP.x2 = x2;
if nargin == 4
  Exact2D = @(t, x, y) DiskExactMode2D(t, x, y, root);
  [~, ~, ~, beta, R0] = Exact2D(0, 0, 0);
  [v1, v2, pr] = L2ProjExact(OP, 0);

  %{
  pf = sqrt(v1.^2 + v2.^2);
  PlotField2D(N, x, y, pf(vmapDG));
  hold on
  PlotMesh2D()
  view(2)
  hold off
  colorbar
  drawnow
  %}

  % Solve Problem
  FinalTime = beta * 2 * pi / R0;

  dt = compute_dt(OP)/2;
  Ntsteps = ceil(FinalTime/dt); dt = FinalTime/Ntsteps;
  [v1,v2,pr,T,err] = Acoustic2D(OP, v1, v2, pr, dt, Ntsteps, 2*p+2);
else
  T = ComputeA(OP);
end
