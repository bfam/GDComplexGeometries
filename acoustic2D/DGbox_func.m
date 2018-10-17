% straight sided simplicial run function
%
% Inputs:
%   p       :: polynomial order is 2*p+1
%   kx      :: (optional; see below) mode number for box solution
%   dt      :: (optional; see below) time step to use
% Outputs:
%   output  :: if nargin == 1 then the matrix is formed and returned
%                 (e.g., dq/dt = A * q)
%              if nargin == 2 then the maximum stable time step for the Taylor
%                 time stepper is determined (looking at energy growth with
%                 random initial condition); Starting time step guess is kx.
%                 output = [dt; dt0] where dt is the found time step and dt0 is
%                 the initial guess
%              if nargin == 3 then the PDE is solved with BoxExactMode2D
%                 solution using kx (dt can be empty in this case and if so then
%                 compute_dt is used to estimate the time step). output the
%                 L2 error at the final time
function [output] = GDbox(p, kx, dt)

addpath ../src
Globals2D_gddg
addpath([NODAL_DG_ROOT, '/Codes1.1/Codes1D'])
addpath([NODAL_DG_ROOT, '/Codes1.1/Codes2D'])
addpath([NODAL_DG_ROOT, '/Codes1.1/ServiceRoutines'])
addpath ../geo/box

Globals2D;

% order of the DG
N = 2*p+1;

% DG initialization
mesh = 'box_r2.msh';
[Nv, VX, VY, K, EToV, BCType] = MeshReaderGmshBC2D(mesh);
StartUp2D;
% BuildBCMaps2D
% build the maps for the boundary conditions and/or mortar
BuildBCMaps2D_gddg()
vmapDG = reshape([1:Np*K], size(x));

% setup the curved grid
intC = 2*(N+1);
OP.dg.cubature = CubatureVolumeMesh2D(intC);
% divide by J^2 because W includes J
OP.dg.cubature.WJI = OP.dg.cubature.W ./ (OP.dg.cubature.J).^2;

intG = 2*(N+1);
OP.dg.gauss = GaussFaceMesh2D(intG);

% build maps for the boundary conditions and mortar with cubature
[OP.dg.gauss.mapB, OP.dg.gauss.mapMor] = ...
  BuildBCMaps2D_cubature(intC, BCType, OP.dg.gauss.mapM);

% Pretend the DG if the last GD block
OP.B{1}.isdg = true;
OP.B{1}.toB = [];
OP.B{1}.toF = [];
OgP.B{1}.quad_order = Nfp-1;
OP.B{1}.vmap = vmapDG;

% Create mortar grids for the DG block

% Here we assume that the "mortar" is aligned with one of the coordinate
% directions with the mapping that the DG faces match the GD faces.

% Set up the global grid and initial conditions
x1 = [x(:)];
x2 = [y(:)];
OP.x1 = x1; OP.x2 = x2;


% Solve Problem
if nargin == 3
  Exact2D = @(t, x, y) BoxExactMode2D(t, x, y, kx);
  FinalTime = 4 / (sqrt(2) * kx);
  if nargin < 3
    dt = [];
  end
  if isempty(dt)
    dt = compute_dt(OP);
  end
  Ntsteps = ceil(FinalTime/dt); dt = FinalTime/Ntsteps;
  [v1, v2, pr] = L2ProjExact(OP, 0);
  [~,~,~,~,err] = Acoustic2D(OP, v1, v2, pr, dt, Ntsteps, 2*p+2);
  output = err;
elseif nargin == 1
  output = ComputeA(OP);
else
  Exact2D = [];
  rng(0);
  v1 = rand(size(x1));
  v2 = rand(size(x1));
  pr = rand(size(x1));

  bracket = [nan, nan];
  dt = compute_dt(OP);
  dt0 = dt;
  if nargin == 2
    dt = kx;
  end

  while sum(isfinite(bracket)) < 2
    [~,~,~,~,~,eng] = Acoustic2D(OP, v1, v2, pr, dt, 100, 2*p+2);
    if eng(end) < eng(1)
      bracket(1) = dt;
      if isfinite(bracket(2))
        break
      else
        dt = 2 * dt;
      end
    else
      bracket(2) = dt;
      if isfinite(bracket(1))
        break
      else
        dt = 0.5 * dt;
      end
    end
  end

  while (bracket(2) - bracket(1)) > 1e-3 * bracket(1)
    dt = (bracket(1) + bracket(2)) / 2;
    [~,~,~,~,~,eng] = Acoustic2D(OP, v1, v2, pr, dt, 100, 2*p+2);
    eng = eng(end) / eng(1);

    if eng < 1
      bracket(1) = dt;
    else
      bracket(2) = dt;
    end
    disp(bracket)
  end

  fprintf('\n\n\nDGbox(p = %d);\n', p);
  fprintf('initial dt    = %e\n', dt0);
  fprintf('final bracket = [%e %e]\n', bracket(1), bracket(2));
  output = [dt; dt0];
end
