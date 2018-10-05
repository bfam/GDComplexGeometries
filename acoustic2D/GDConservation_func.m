% Conservation test code
% Inputs:
%   p       :: GD polynomial order is 2*p+1
%   N0      :: GD grid uses N0+1 grid points
%   finterp :: 0 = use L2 projected coordinates
%           :: 1 = make curved faces coordinates to be polynomial order 2*p+1
%           ::     after L2 projecting the coordinates
%           :: 2 = resample the curved faces after L2 projecting the
%           ::     coordinates (this ensures that the conforming curved meshes
%           ::     have same surface metric terms)
%   Ng      :: Number of free ghost points to use
%   alpha   :: (optional) upwinding parameter: 1 -> upwind (default)
%           ::                                 0 -> central
% Outputs:
%   T    :: time
%   eng  :: energy
%   m_pr :: L2 norm of pressure
%   m_v1 :: L2 norm of v1
%   m_v2 :: L2 norm of v2
function [T, eng, m_pr, m_v1, m_v2] = GDConservation_func(p, N0, finterp, ...
                                                          Ng, alpha)

if nargin < 5
  alpha = 1;
end

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
if finterp == 1
  grid{1}.interp = [ 0, 1,  0, 1];
elseif finterp == 2
  grid{1}.interp = [ 0, 0,  0, 0];
end

grid{2}.x1 = @(r1, r2) x((r1+1)/2, (r2-1)/2);
grid{2}.x2 = @(r1, r2) y((r1+1)/2, (r2-1)/2);
grid{2}.N1 = N0(2); grid{2}.N2 = N0(2);
if finterp == 1
  grid{2}.interp = [1,  0,  0, 1];
elseif finterp == 2
  grid{2}.interp = [ 0, 0,  0, 0];
end

grid{3}.x1 = @(r1, r2) x((r1-1)/2, (r2+1)/2);
grid{3}.x2 = @(r1, r2) y((r1-1)/2, (r2+1)/2);
grid{3}.N1 = N0(3); grid{3}.N2 = N0(3);
if finterp == 1
  grid{3}.interp = [ 0, 1, 1,  0];
elseif finterp == 2
  grid{3}.interp = [ 0, 0,  0, 0];
end

grid{4}.x1 = @(r1, r2) x((r1+1)/2, (r2+1)/2);
grid{4}.x2 = @(r1, r2) y((r1+1)/2, (r2+1)/2);
grid{4}.N1 = N0(4); grid{4}.N2 = N0(4);
if finterp == 1
  grid{4}.interp = [1,  0, 1,  0];
elseif finterp == 2
  grid{4}.interp = [ 0, 0,  0, 0];
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
if alpha == 1
  dt = compute_dt(OP) / 2;
else
  dt = compute_dt(OP) / 10;
end

% clear this so we do not try to calate error
clear Exact2D

% random initial condition
rng(1);
pr = rand(size(x1));
v1 = rand(size(x1));
v2 = rand(size(x1));

gd_contourf(OP.B, pr, true, 100);colorbar
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
OP.alpha = alpha;
[~,~,~,T,~,eng,~,~,m_pr,m_v1,m_v2] = ...
  Acoustic2D(OP, v1, v2, pr, dt, Ntsteps, 2*abs(p)+2, true);
