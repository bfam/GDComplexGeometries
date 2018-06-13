% Curved, noncoforming GD Box test
% Inputs:
%   p       :: GD polynomial order is 2*p+1
%   N0      :: GD grid uses N0+1 grid points
%   kx      :: mode number for box solution
%   Ng      :: (optional) Number of free ghost points to use (default: p)
% Outputs:
%   T    :: time
%   err  :: energy
function [T, err] = GDbox_curved_func(p, N0, kx, Ng)

addpath ../src
Globals2D_gddg;

if length(N0) == 1
  N0 = [N0, N0, N0, N0];
end
if nargin < 4
  Ng = [];
end
if isempty(Ng)
  Ng = p;
end


% Set up the grid structure for the four blocks
% +---+---+
% | 3 | 4 |
% +---+---+
% | 1 | 2 |
% +---+---+
beta_in = pi / 4;
beta = @(r,s) beta_in * (1-r.^2) .* (1-s.^2);
x = @(r,s)  r .* cos(beta(r,s)) - s .* sin(beta(r,s));
y = @(r,s)  r .* sin(beta(r,s)) + s .* cos(beta(r,s));

grid{1}.x1 = @(r1, r2) x((r1-1)/2, (r2-1)/2);
grid{1}.x2 = @(r1, r2) y((r1-1)/2, (r2-1)/2);
grid{1}.lx = [-1 0]; grid{1}.ly = [-1 0];
grid{1}.N1 = N0(1); grid{1}.N2 = N0(1);

grid{2}.x1 = @(r1, r2) x((r1+1)/2, (r2-1)/2);
grid{2}.x2 = @(r1, r2) y((r1+1)/2, (r2-1)/2);
grid{2}.lx = [0 1]; grid{2}.ly = [-1 0];
grid{2}.N1 = N0(2); grid{2}.N2 = N0(2);

grid{3}.x1 = @(r1, r2) x((r1-1)/2, (r2+1)/2);
grid{3}.x2 = @(r1, r2) y((r1-1)/2, (r2+1)/2);
grid{3}.lx = [-1 0]; grid{3}.ly = [0 1];
grid{3}.N1 = N0(3); grid{3}.N2 = N0(3);

grid{4}.x1 = @(r1, r2) x((r1+1)/2, (r2+1)/2);
grid{4}.x2 = @(r1, r2) y((r1+1)/2, (r2+1)/2);
grid{4}.lx = [0 1]; grid{4}.ly = [0 1];
grid{4}.N1 = N0(4); grid{4}.N2 = N0(4);

% Create the GD element operators
st = 0; x1 = []; x2 = [];
for b = 1:length(grid)
  OP.B{b} = gd_setup_curved_alias(grid{b}.N1, grid{b}.N2, p, 2*(p+2), grid{b}, Ng);

  % Set up the global vmap
  OP.B{b}.vmap = st + (1:length(OP.B{b}.x1))';
  st = OP.B{b}.vmap(end);
  OP.B{b}.vmap = reshape(OP.B{b}.vmap, OP.B{b}.Np2, OP.B{b}.Np1);

  % Set up the global grid
  % plot(OP.B{b}.x1, OP.B{b}.x2, '*')
  % hold on
  x1 = [x1;OP.B{b}.x1]; x2 = [x2;OP.B{b}.x2];
end
OP.x1 = x1; OP.x2 = x2;
% hold off
% axis image

% set up the coupling
% zeros means boundary conditions
% negative face means flipped face
OP.B{1}.toB = [0 2 0 3];
OP.B{1}.toF = [2 1 4 3];

OP.B{2}.toB = [1 0 0 4];
OP.B{2}.toF = [2 1 4 3];

OP.B{3}.toB = [0 4 1 0];
OP.B{3}.toF = [2 1 4 3];

OP.B{4}.toB = [3 0 2 0];
OP.B{4}.toF = [2 1 4 3];

OP.B = mortar_create(OP.B);
% plot_mesh(OP.B);

% Estimate time step
dt = compute_dt(OP)/2;

% initial condition is the (approximate) L2 projection of the exact solution
Exact2D = @(t, x, y) BoxExactMode2D(t, x, y, kx);
[v1, v2, pr] = L2ProjExact(OP, 0);

% run the simulation
FinalTime = 4 / (sqrt(2) * kx);
Ntsteps = ceil(FinalTime/dt); dt = FinalTime/Ntsteps;
[v1,v2,pr,T,err] = Acoustic2D(OP, v1, v2, pr, dt, Ntsteps, 2*p+2);
