% straight sided GD run function
%
% Inputs:
%   p       :: polynomial order is 2*p+1
%   Ng      :: number of free ghost points
%   N0      :: N0+1 interior GD grid points (in each dimension)
%   kx      :: (optional; see below) mode number for box solution
%   dt      :: (optional; see below) time step to use
% Outputs:
%   output  :: if nargin == 3 then the matrix is formed and returned
%                 (e.g., dq/dt = A * q)
%              if nargin == 4 then the maximum stable time step for the Taylor
%                 time stepper is determined (looking at energy growth with
%                 random initial condition); Starting time step guess is kx.
%                 output = [dt; dt0] where dt is the found time step and dt0 is
%                 the initial guess
%              if nargin == 5 then the PDE is solved with BoxExactMode2D
%                 solution using kx (dt can be empty in this case and if so then
%                 compute_dt is used to estimate the time step). output the
%                 L2 error at the final time
function output = GDbox(p, Ng, N0, kx, dt)

addpath ../src
addpath ../geo/box

Globals2D_gddg

% order of gd approximation

% order of the mortar
Nm = 2*p+1;

%
% Set up the Geometry transforms
%

% center block
st = 0; x1 = []; x2 = [];
OP.B{1} = gd_setup_affine(N0, N0, p, Nm+1, [-1, 1], [-1 1], Ng);

% Set up the global vmap
OP.B{1}.vmap = st + (1:length(OP.B{1}.x1))';
st = OP.B{1}.vmap(end);
OP.B{1}.vmap = reshape(OP.B{1}.vmap, OP.B{1}.Np2, OP.B{1}.Np1);

% Set up the global grid
x1 = [x1;OP.B{1}.x1]; x2 = [x2;OP.B{1}.x2];
OP.x1 = x1; OP.x2 = x2;


OP.B{1}.toB = [0 0 0 0];
OP.B{1}.toF = [1 2 3 4];

% Solve Problem
if nargin == 3
  output = ComputeA(OP);
elseif nargin == 5
  Exact2D = @(t, x, y) BoxExactMode2D(t, x, y, kx);
  FinalTime = 4 / (sqrt(2) * kx);
  if nargin < 5
    dt = [];
  end
  if isempty(dt)
    dt = compute_dt(OP);
  end
  Ntsteps = ceil(FinalTime/dt); dt = FinalTime/Ntsteps;
  [v1, v2, pr] = L2ProjExact(OP, 0);
  [~,~,~,~,err] = Acoustic2D(OP, v1, v2, pr, dt, Ntsteps, 2*p+2);
  output = err;
else
  Exact2D = [];
  rng(0);
  v1 = rand(size(x1));
  v2 = rand(size(x1));
  pr = rand(size(x1));

  bracket = [nan, nan];
  dt = compute_dt(OP);
  dt0 = dt;
  if nargin == 4
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

  fprintf('\n\n\nGDbox(p = %d, Ng = %d, N0 = %d);\n', p, Ng, N0);
  fprintf('initial dt    = %e\n', dt0);
  fprintf('final bracket = [%e %e]\n', bracket(1), bracket(2));
  output = [dt; dt0];
end
