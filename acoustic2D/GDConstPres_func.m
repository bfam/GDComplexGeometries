% Conservation test code
% Inputs:
%   p       :: GD polynomial order is 2*p+1
%   N0      :: GD grid uses N0+1 grid points
%   do_proj :: 0 = use L2 projected coordinates
%           :: 1 = make curved faces coordinates to be polynomial order 2*p+1
%           ::     after L2 projecting the coordinates
%           :: 2 = resample the curved faces after after L2 projecting the
%           ::     coordinates (this ensures that the conforming curved meshes
%           ::     have same surface metric terms)
%   Ng      :: Number of free ghost points to use
% Outputs:
%   T    :: time
%   err  :: error
function [T, err] = GDConstPres_func(p, N0, do_proj, Ng)

addpath ../src
Globals2D_gddg;

% Global coordinate transform
beta_in = pi / 4;
beta = @(r,s) beta_in * (1-r.^2) .* (1-s.^2);
x = @(r,s)  r .* cos(beta(r,s)) + s .* sin(beta(r,s));
y = @(r,s) -r .* sin(beta(r,s)) + s .* cos(beta(r,s));


% Set up the grid structure for the four blocks
% +---+---+
% | 3 | 4 |
% +---+---+
% | 1 | 2 |
% +---+---+
grid{1}.x1 = @(r1, r2) x((r1-1)/2, (r2-1)/2);
grid{1}.x2 = @(r1, r2) y((r1-1)/2, (r2-1)/2);
grid{1}.N1 = N0(1); grid{1}.N2 = N0(1);
if do_proj == 1
  grid{1}.proj = [ 0, 20,  0, 20];
elseif do_proj == 2
  grid{1}.proj = [ 0, 0,  0, 0];
end

grid{2}.x1 = @(r1, r2) x((r1+1)/2, (r2-1)/2);
grid{2}.x2 = @(r1, r2) y((r1+1)/2, (r2-1)/2);
grid{2}.N1 = N0(2); grid{2}.N2 = N0(2);
if do_proj == 1
  grid{2}.proj = [20,  0,  0, 20];
elseif do_proj == 2
  grid{2}.proj = [ 0, 0,  0, 0];
end

grid{3}.x1 = @(r1, r2) x((r1-1)/2, (r2+1)/2);
grid{3}.x2 = @(r1, r2) y((r1-1)/2, (r2+1)/2);
grid{3}.N1 = N0(3); grid{3}.N2 = N0(3);
if do_proj == 1
  grid{3}.proj = [ 0, 20, 20,  0];
elseif do_proj == 2
  grid{3}.proj = [ 0, 0,  0, 0];
end

grid{4}.x1 = @(r1, r2) x((r1+1)/2, (r2+1)/2);
grid{4}.x2 = @(r1, r2) y((r1+1)/2, (r2+1)/2);
grid{4}.N1 = N0(4); grid{4}.N2 = N0(4);
if do_proj == 1
  grid{4}.proj = [20,  0, 20,  0];
elseif do_proj == 2
  grid{4}.proj = [ 0, 0,  0, 0];
end

% Create the GD element operators
st = 0; x1 = []; x2 = [];
for k = 1:length(grid)
  OP.B{k} = gd_setup_curved_alias(grid{k}.N1, grid{k}.N2, p, 2*(p+2), grid{k}, Ng);

  % Set up the global vmap
  OP.B{k}.vmap = st + (1:length(OP.B{k}.x1))';
  st = OP.B{k}.vmap(end);
  OP.B{k}.vmap = reshape(OP.B{k}.vmap, OP.B{k}.Np2, OP.B{k}.Np1);

  % Set up the global grid
  plot(OP.B{k}.x1, OP.B{k}.x2, '*')
  hold on
  x1 = [x1;OP.B{k}.x1]; x2 = [x2;OP.B{k}.x2];
end
OP.x1 = x1; OP.x2 = x2;
hold off
axis image

% set up the coupling
% zeros means boundary conditions
% negative face means flipped face
OP.B{1}.toB = [2 2 3 3];
OP.B{1}.toF = [2 1 4 3];
OP.B{1}.toShftX = [2 0 0 0];
OP.B{1}.toShftY = [0 0 2 0];

OP.B{2}.toB = [1 1 4 4];
OP.B{2}.toF = [2 1 4 3];
OP.B{2}.toShftX = [0 -2 0 0];
OP.B{2}.toShftY = [0  0 2 0];

OP.B{3}.toB = [4 4 1 1];
OP.B{3}.toF = [2 1 4 3];
OP.B{3}.toShftX = [2 0 0  0];
OP.B{3}.toShftY = [0 0 0 -2];

OP.B{4}.toB = [3 3 2 2];
OP.B{4}.toF = [2 1 4 3];
OP.B{4}.toShftX = [0 -2 0  0];
OP.B{4}.toShftY = [0  0 0 -2];

OP.B = mortar_create(OP.B);
plot_mesh(OP.B);

% Estimate time step
dt = compute_dt(OP)/2;

% Set up constant initial conditions
Exact2D = @(t, x, y) deal(1 * ones(size(x)), 2 * ones(size(x)), 3 * ones(size(x)));
[v1, v2, pr0] = Exact2D(0, x1, x2);
%{
hold on
for k = 1:length(grid)
  OP.B{k}.x1 = OP.B{k}.x1 + 2;
end
gd_contourf(OP.B, pr0, true, 100);colorbar
hold off
return
hold off
%}
drawnow

% Run the simulation
FinalTime = 1;
Ntsteps = ceil(FinalTime/dt); dt = FinalTime/Ntsteps;
[v1,v2,pr,T,err] = Acoustic2D(OP, v1, v2, pr0, dt, Ntsteps, 2*p+2);
gd_contourf(OP.B, pr-pr0, true, 100);colorbar
