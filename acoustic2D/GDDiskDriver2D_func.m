% Curved, noncoforming GD Disk test
% Inputs:
%   p       :: GD polynomial order is 2*p+1
%   N0      :: GD grid uses N0+1 grid points
%   sr      :: level of refiment, e.g, N-> N0*2^sr
%   root    :: (optional) which root for the solution to use (default: [])
%              NOTE: if isempty(root) then the codes does not solve the disk
%                    problem but instead computes the maximum stable timestep
%                    using bisection on the time step size and a random initial
%                    solution
%   Ng      :: (optional) Number of free ghost points to use (default: p)
% Outputs:
%   T    :: time
%   err  :: energy
function [T, err] = GDDiskDriver2D_func(p, N0, sr, root, Ng)

if nargin < 4
  root = [];
end
if nargin < 5
  Ng = [];
end
if isempty(Ng)
  Ng = p;
end

if isempty(Ng)
  fprintf('+--------------------------------------------------------------+\n');
  fprintf('| GDDiskDriver2D_func(p = %3d, N0 = %3d, sr = %3d, root = %3d) |\n', p, N0, sr, root)
  fprintf('+--------------------------------------------------------------+\n');
else
  fprintf('+----------------------------------------------------------------------+\n');
  fprintf('| GDDiskDriver2D_func(p = %3d, N0 = %3d, sr = %3d, root = %3d Ng = %3d) |\n', ...
          p, N0, sr, Ng, root)
  fprintf('+----------------------------------------------------------------------+\n');
end

addpath ../src
addpath ../geo/disk

Globals2D_gddg;

% order of the mortar
Nm = 2*p+1;

%
% Set up the Geometry transforms
%

% center block
grid{1}.x1 = @(r1, r2) disktransform(r1, r2, 1, 1);
grid{1}.x2 = @(r1, r2) disktransform(r1, r2, 1, 2);
grid{1}.N1 = N0; grid{1}.N2 = N0;

h0 = (grid{1}.x1(1,0) - grid{1}.x1(-1,0)) / N0;
N1 = N0;
h1 = max(diff(disktransform(linspace(-1,1,N1+1), 0, 2, 1)));
disp([h0, h1])
while h1 < h0
  N1 = N1-1;
  h1 = max(diff(disktransform(linspace(-1,1,N1+1), 0, 2, 1)));
end
while h1 > h0
  N1 = N1+1;
  h1 = max(diff(disktransform(linspace(-1,1,N1+1), 0, 2, 1)));
end
N2 = ceil((pi/2)/h0);


% left block
grid{2}.x1 = @(r1, r2) disktransform(r1, r2, 2, 1);
grid{2}.x2 = @(r1, r2) disktransform(r1, r2, 2, 2);
grid{2}.N1 = N1; grid{2}.N2 = N2;

% right block
grid{3}.x1 = @(r1, r2) disktransform(r1, r2, 3, 1);
grid{3}.x2 = @(r1, r2) disktransform(r1, r2, 3, 2);
grid{3}.N1 = N1; grid{3}.N2 = N2;

% bottom block
grid{4}.x1 = @(r1, r2) disktransform(r1, r2, 4, 1);
grid{4}.x2 = @(r1, r2) disktransform(r1, r2, 4, 2);
grid{4}.N1 = N2; grid{4}.N2 = N1;

% top block
grid{5}.x1 = @(r1, r2) disktransform(r1, r2, 5, 1);
grid{5}.x2 = @(r1, r2) disktransform(r1, r2, 5, 2);
grid{5}.N1 = N2; grid{5}.N2 = N1;


for k = 1:length(grid)
  grid{k}.N1 = grid{k}.N1 * 2^sr;
  grid{k}.N2 = grid{k}.N2 * 2^sr;
  fprintf('grid{%d} (Nx, Ny) = (%d, %d)\n', k, grid{k}.N1, grid{k}.N2);
end

%
% Initialize the blocks
%

st = 0; x1 = []; x2 = [];
for k = 1:length(grid)
  if k == 1
    OP.B{k} = gd_setup_affine(grid{k}.N1, grid{k}.N2, p, Nm+1, ...
                              [-1 1]/3, [-1 1]/3, Ng);
  else
    OP.B{k} = gd_setup_curved_alias(grid{k}.N1, grid{k}.N2, p, Nm+1, grid{k}, Ng);
  end

  % Set up the global vmap
  OP.B{k}.vmap = st + (1:length(OP.B{k}.x1))';
  st = OP.B{k}.vmap(end);
  OP.B{k}.vmap = reshape(OP.B{k}.vmap, OP.B{k}.Np2, OP.B{k}.Np1);

  % Set up the global grid
  x1 = [x1;OP.B{k}.x1]; x2 = [x2;OP.B{k}.x2];
end
OP.x1 = x1; OP.x2 = x2;

% set up the coupling
% zeros means boundary conditions
% negative face means flipped face
OP.B{1}.toB = [2 3 4 5];
OP.B{1}.toF = [2 1 4 3];

OP.B{2}.toB = [0 1 4 5];
OP.B{2}.toF = [1 1 1 1];

OP.B{3}.toB = [1 0 4 5];
OP.B{3}.toF = [2 2 2 2];

OP.B{4}.toB = [2 3 0 1];
OP.B{4}.toF = [3 3 3 3];

OP.B{5}.toB = [2 3 1 0];
OP.B{5}.toF = [4 4 4 4];

OP.B = mortar_create(OP.B);

if ~isempty(root) % Solve PDE
  Exact2D = @(t, x, y) DiskExactMode2D(t, x, y, root);
  [~, ~, ~, beta, R0] = Exact2D(0, 0, 0);
  [v1, v2, pr] = L2ProjExact(OP, 0);

  % gd_contourf(OP.B, pr, false, 100);
  % gd_contourf(OP.B, sqrt(v1.^2+v2.^2), false, 100);
  % drawnow

  % Solve Problem
  FinalTime = beta * 2 * pi / R0;

  dt = compute_dt(OP)/2;
  Ntsteps = ceil(FinalTime/dt); dt = FinalTime/Ntsteps;
  [v1,v2,pr,T,err] = Acoustic2D(OP, v1, v2, pr, dt, Ntsteps, 2*p+2);
  % A = AcousticCreateMatrix(OP);
  % spy(A);
else % Compute maximum time step
  dt = compute_dt(OP)/2;
  Exact2D = [];
  rng(0);
  v1 = rand(size(x1));
  v2 = rand(size(x1));
  pr = rand(size(x1));

  bracket = [nan, nan];
  dt = compute_dt(OP);
  dt0 = dt;

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
  T = dt;
end
