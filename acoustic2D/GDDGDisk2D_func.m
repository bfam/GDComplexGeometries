% Disk test problem with couple affine GD elements and curved simplicial
% elements
%
% Inputs:
%   p         :: GD polynomial order is 2*p+1
%   N0        :: GD grid uses N0+1 grid points
%   sr        :: level of refiment, e.g, N-> N0*2^sr
%   mesh_base :: mesh used is sprintf('%s_r%d.msh', mesh_base, sr)
%   root      :: (optional) which root for the solution to use (default: [])
%                NOTE: if isempty(root) then the codes does not solve the disk
%                      problem but instead computes the operator matrix
%   Ng       :: (optional) Number of free ghost points to use (default: p)
% Outputs:
%   T    :: time
%   err  :: energy
function [T,err] = GDDGDisk2D_func(p, N0, sr, mesh_base, root, Ng)
if nargin < 5
  root = [];
end
if nargin < 6
  Ng = [];
end
if isempty(Ng)
  Ng = p;
end

fprintf('GDDGDisk2D_func(%3d, %3d, %3d, %s)\n', p, N0, sr, mesh_base)

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
intC = 2*(N+1);
OP.dg.cubature = DGMetricStorageUpdate(CubatureVolumeMesh2D(intC));
%{
OP.dg.cubature = CubatureVolumeMesh2D(intC);
OP.dg.cubature.WJI = OP.dg.cubature.W ./ (OP.dg.cubature.J).^2;
%}

intG = 2*(N+1);
OP.dg.gauss = GaussFaceMesh2D(intG);
[OP.dg.gauss.mapB, OP.dg.gauss.mapMor] = ...
  BuildBCMaps2D_cubature(intG, BCType, OP.dg.gauss.mapM);

% GD initialization
OP.B{1} = gd_setup_affine(N0 * 2^sr, N0 * 2^sr, p, Nm+1, [-1 1]*sqrt(2)/2, [-1 1]*sqrt(2)/2, Ng);

% Set up the global vmap
OP.B{1}.vmap = vmapDG(end) + (1:length(OP.B{1}.x1))';
if isfield(OP.B{1}, 'isaffine')
  OP.B{1}.vmap = reshape(OP.B{1}.vmap, OP.B{1}.Np2, OP.B{1}.Np1);
end

OP.B{1}.toB = [2 2 2 2];
OP.B{1}.toF = [1 2 3 4];

% Pretend the DG if the last GD block
OP.B{2}.isdg = true;
OP.B{2}.toB = [1 1 1 1];
OP.B{2}.toF = [1 2 3 4];
OP.B{2}.quad_order = Nfp-1;
OP.B{2}.vmap = vmapDG;

% Create mortar grids for the DG block

% Here we assume that the "mortar" is aligned with one of the coordinate
% directions with the mapping that the DG faces match the GD faces.
for m = 1:length(vmapMor)
  clear f

  % First we sort the maps so that they match the expected GD ordering
  for k = 1:size(vmapMor{m},2)
    if m == 1 || m == 2
      [~, I] = sort(y(vmapMor{m}(:,k)));
    else
      [~, I] = sort(x(vmapMor{m}(:,k)));
    end
    vmapMor{m}(:,k) = vmapMor{m}(I,k);
    mapMor{m}(:,k) = mapMor{m}(I,k);
  end
  if m == 1 || m == 2
     [~, I] = sort(y(vmapMor{m}(1,:)));
  else
     [~, I] = sort(x(vmapMor{m}(1,:)));
  end
  vmapMor{m} = vmapMor{m}(:,I);
  mapMor{m} = mapMor{m}(:,I);

  % Now we create the data needed to create the mortar
  % (see gd_setup_curved)
  f.x1 = x(vmapMor{m}(:));
  f.x2 = y(vmapMor{m}(:));
  f.corners = [f.x1(  1,   1), f.x2(  1,   1);
               f.x1(end, end), f.x2(end, end)];
  if m == 1 || m == 2
     f.g  = [y(vmapMor{m}(1,:)),y(vmapMor{m}(end,end))]';
     f.rq = f.x2(:);
  else
     f.g = [x(vmapMor{m}(1,:)),x(vmapMor{m}(end,end))]';
     f.rq = f.x1(:);
  end
  f.rq = 2*(f.rq-f.g(1))/(f.g(end)-f.g(1))-1;
  f.g  = 2*(f.g -f.g(1))/(f.g(end)-f.g(1))-1;
  f.P = 1;
  f.vmap = vmapMor{m}(:);

  OP.B{2}.f{m} = f;
end

OP.B = mortar_create(OP.B);

% Set up the global grid and initial conditions
x1 = [x(:);OP.B{1}.x1];
x2 = [y(:);OP.B{1}.x2];
OP.x1 = x1; OP.x2 = x2;

if ~isempty(root)
  Exact2D = @(t, x, y) DiskExactMode2D(t, x, y, root);
  [~, ~, ~, beta, R0] = Exact2D(0, 0, 0);
  [v1, v2, pr] = L2ProjExact(OP, 0);

  %{
  pf = sqrt(v1.^2 + v2.^2);
  gd_contourf(OP.B, pf, true, 100);
  hold on
  PlotField2D(N, x, y, pf(vmapDG));
  PlotMesh2D()
  hold off
  colorbar
  drawnow
  return
  %}

  % Solve Problem
  FinalTime = beta * 2 * pi / R0;

  dt = compute_dt(OP)/2;
  Ntsteps = ceil(FinalTime/dt); dt = FinalTime/Ntsteps;
  [v1,v2,pr,T,err] = Acoustic2D(OP, v1, v2, pr, dt, Ntsteps, 2*p+2);

  %{
  pf = sqrt(v1.^2 + v2.^2);
  gd_contourf(OP.B, pf, true, 100);
  hold on
  PlotField2D(N, x, y, pf(vmapDG));
  PlotMesh2D()
  hold off
  colorbar
  drawnow
  %}
else
  T = ComputeA(OP);
end
