% Setup curved GD operator
% Inputs:
%   N1             :: Interior grid in the r1 direction is size N1+1
%   N2             :: Interior grid in the r2 direction is size N2+1
%   p              :: GD polynomial order is 2*p+1
%   quad_order     :: order of quadrature to use for subcells
%   lim1, lim2     :: limits for x1 and x2
%   Ng             :: (optional) Number of free ghost points [default Ng = p]
% Outputs:
%   B :: GD operator struct
% members of the struct:
%   isaffine                :: boolean saying that this GD element is affine
%   N1, N2                  :: interior grid dimensions
%   p                       :: polynomial order is 2*p+1
%   Ng                      :: number of ghost points
%   Np1, Np2                :: total number of points in each dimension
%   quad_order              :: subcell quadrature order
%   P1_ge, P2_ge            :: projections to equally spaced grid
%   M1_1d, M2_1d            :: 1-D mass matrices
%   S1_1d, S2_1d            :: 1-D stiffness matrices
%   J                       :: element Jacobian determinant (constant)
%   P1, P2                  :: projections to equally spaced grid
%   w1, w2                  :: 1D quadrature grids
%   w                       :: 2D quadrature grid
%   E1, E2                  :: non-free ghost extrapolations
%   P2F, P1F                :: 2D quadrature interpolation [apply: P1F*(P2F*q)]
%   h                       :: grid spacing in physical coordinates
%   x1, x2                  :: 2D GD coordinate points
%   xq1, xq2                :: physicaly coordinates of quadrature points
%   M                       :: 2D GD mass matrix
%   R1_1d, R2_1d            :: 1D mass matrix Cholesky factors
%   Sx, Sy                  :: functions to apply stiffness matrices
%   SxT, SyT                :: functions to apply stiffness matrices transpose
%   massInv                 :: function to apply mass matrix inverse
%   f                       :: face map struct
%     -> vmap               :: volume points on this face
%     -> P                  :: interpolation to face quadrature
%     -> w                  :: face quadrature weights
%     -> g                  :: face cell edges
%     -> rq                 :: face quadrature points
%     -> x1, x2             :: physical coordinates of face GD points
%     -> sJ                 :: surface Jacobian
%     -> n1, n2             :: surface unit normal
function [B] = gd_setup_affine(N1, N2, p, quad_order, lim1, lim2, Ng, err_quad_order)

  if nargin < 7
    Ng = [];
  end
  if isempty(Ng)
    Ng = p;
  end
  if nargin < 8
    err_quad_order = [];
  end
  if isempty(err_quad_order)
    err_quad_order = quad_order;
  end

  N_plot_points = 10;

  B.isaffine = true;

  B.N1  = N1;
  B.N2  = N2;
  B.p   = p;
  B.Ng  = Ng;

  % number of points with ghost
  Ng1 = B.N1 + 1 + 2 * B.p;
  Ng2 = B.N2 + 1 + 2 * B.p;
  B.Np1 = B.N1 + 1 + 2 * Ng;
  B.Np2 = B.N2 + 1 + 2 * Ng;

  % compute the gd quadrature rule
  B.quad_order = quad_order;
  [P1_1d, w1_1d, r1_1d, rq1_1d, D1_1d, B.P1_ge] = ...
    gd_quadrature(p, B.quad_order, N1, Ng1, N_plot_points);
  sc1 = lim1(2) - lim1(1);
  B.M1_1d = (sc1 / 2) * P1_1d' * diag(sparse(w1_1d)) * P1_1d;
  B.S1_1d = P1_1d' * diag(sparse(w1_1d)) * D1_1d * P1_1d;

  [P2_1d, w2_1d, r2_1d, rq2_1d, D2_1d, B.P2_ge] = ...
    gd_quadrature(p, B.quad_order, N2, Ng2, N_plot_points);
  sc2 = lim2(2) - lim2(1);
  B.M2_1d = (sc2 / 2) * P2_1d' * diag(sparse(w2_1d)) * P2_1d;
  B.S2_1d = P2_1d' * diag(sparse(w2_1d)) * D2_1d * P2_1d;

  B.J = sc2 * sc1 / 4;

  B.P1 = P1_1d;
  B.P2 = P2_1d;


  B.w1 = w1_1d;
  B.w2 = w2_1d;
  B.w  = kron(B.w1, B.w2);

  E1 = ghost_extrapolation(N1, Ng, p);
  E2 = ghost_extrapolation(N2, Ng, p);
  B.E1 = E1;
  B.E2 = E2;
  B.P2F = kron(speye(size(B.P1*E1,2)), B.P2*E2);
  B.P1F = kron(B.P1*E1, speye(size(B.P2*E2,1)));

  B.P1_ge = B.P1_ge * E1;
  B.P2_ge = B.P2_ge * E2;

  B.M1_1d = E1' * B.M1_1d * E1;
  B.S1_1d = E1' * B.S1_1d * E1;
  B.M2_1d = E2' * B.M2_1d * E2;
  B.S2_1d = E2' * B.S2_1d * E2;

  P1_1d = P1_1d * E1;
  P2_1d = P2_1d * E2;

  B.h = min(sc1/N1,sc2/N2);

  % Set up reference grid
  x1 = lim1(2) * (1 + r1_1d) / 2 + lim1(1) * (1 - r1_1d) / 2;
  x2 = lim2(2) * (1 + r2_1d) / 2 + lim2(1) * (1 - r2_1d) / 2;
  x1 = x1(p-Ng+1:end-(p-Ng));
  x2 = x2(p-Ng+1:end-(p-Ng));
  [x1,x2] = meshgrid(x1, x2);

  xq1 = lim1(2) * (1 + rq1_1d) / 2 + lim1(1) * (1 - rq1_1d) / 2;
  xq2 = lim2(2) * (1 + rq2_1d) / 2 + lim2(1) * (1 - rq2_1d) / 2;
  [xq1,xq2] = meshgrid(xq1, xq2);
  B.xq1 = xq1;
  B.xq2 = xq2;

  B.x1 = x1(:);
  B.x2 = x2(:);

  B.M  = kron(B.M1_1d, B.M2_1d);
  B.R1_1d  = chol(B.M1_1d);
  B.R2_1d  = chol(B.M2_1d);


  B.Sx  = @(v) B.M2_1d  * (v * B.S1_1d');
  B.Sy  = @(v) B.S2_1d  * (v * B.M1_1d );
  B.SxT = @(v) B.M2_1d  * (v * B.S1_1d );
  B.SyT = @(v) B.S2_1d' * (v * B.M1_1d );

  B.massInv = @(v) (B.R2_1d \ (B.R2_1d' \ ((v / B.R1_1d) / B.R1_1d')));

  %
  % Face stuff
  %

  % Boundary maps
  B.f{1}.vmap  =  Ng      *B.Np2 +       (1:B.Np2)';
  B.f{2}.vmap  = (Ng+B.N1)*B.Np2 +       (1:B.Np2)';
  B.f{3}.vmap = Ng+1           + B.Np2*(0:B.Np1-1)';
  B.f{4}.vmap = Ng+B.N2+1      + B.Np2*(0:B.Np1-1)';

  B.f{1}.P = P2_1d; B.f{1}.w = w2_1d;
  B.f{2}.P = P2_1d; B.f{2}.w = w2_1d;
  B.f{3}.P = P1_1d; B.f{3}.w = w1_1d;
  B.f{4}.P = P1_1d; B.f{4}.w = w1_1d;

  %
  % Boundary Metric Terms
  %
  B.f{1}.g  = linspace(-1,1,B.N2+1)';                % gd cells edges
  B.f{1}.x1 = lim1(1) * ones(size(r2_1d));
  B.f{1}.x2 = lim2(2) * (1 + r2_1d) / 2 + lim2(1) * (1 - r2_1d) / 2;
  B.f{1}.sJ = (lim2(2) - lim2(1)) / 2;
  B.f{1}.n1 = -1;
  B.f{1}.n2 =  0;
  B.f{1}.x1 = B.f{1}.x1(p-Ng+1:end-(p-Ng));
  B.f{1}.x2 = B.f{1}.x2(p-Ng+1:end-(p-Ng));

  B.f{2}.g  = linspace(-1,1,B.N2+1)';                % gd cells edges
  B.f{2}.x1 = lim1(2) * ones(size(r2_1d));
  B.f{2}.x2 = lim2(2) * (1 + r2_1d) / 2 + lim2(1) * (1 - r2_1d) / 2;
  B.f{2}.sJ = (lim2(2) - lim2(1)) / 2;
  B.f{2}.n1 = 1;
  B.f{2}.n2 = 0;
  B.f{2}.x1 = B.f{2}.x1(p-Ng+1:end-(p-Ng));
  B.f{2}.x2 = B.f{2}.x2(p-Ng+1:end-(p-Ng));

  B.f{3}.g  = linspace(-1,1,B.N1+1)';                % gd cells edges
  B.f{3}.x1 = lim1(2) * (1 + r1_1d) / 2 + lim1(1) * (1 - r1_1d) / 2;
  B.f{3}.x2 = lim2(1) * ones(size(r1_1d));
  B.f{3}.sJ = (lim1(2) - lim1(1)) / 2;
  B.f{3}.n1 =  0;
  B.f{3}.n2 = -1;
  B.f{3}.x1 = B.f{3}.x1(p-Ng+1:end-(p-Ng));
  B.f{3}.x2 = B.f{3}.x2(p-Ng+1:end-(p-Ng));

  B.f{4}.g  = linspace(-1,1,B.N1+1)';                % gd cells edges
  B.f{4}.x1 = lim1(2) * (1 + r1_1d) / 2 + lim1(1) * (1 - r1_1d) / 2;
  B.f{4}.x2 = lim2(2) * ones(size(r1_1d));
  B.f{4}.sJ = (lim1(2) - lim1(1)) / 2;
  B.f{4}.n1 = 0;
  B.f{4}.n2 = 1;
  B.f{4}.x1 = B.f{4}.x1(p-Ng+1:end-(p-Ng));
  B.f{4}.x2 = B.f{4}.x2(p-Ng+1:end-(p-Ng));

  B.f{1}.corners = [lim1(1), lim2(1);lim1(1), lim2(2)];
  B.f{2}.corners = [lim1(2), lim2(1);lim1(2), lim2(2)];
  B.f{3}.corners = [lim1(1), lim2(1);lim1(2), lim2(1)];
  B.f{4}.corners = [lim1(1), lim2(2);lim1(2), lim2(2)];

  B.f{1}.rq = rq2_1d;
  B.f{2}.rq = rq2_1d;
  B.f{3}.rq = rq1_1d;
  B.f{4}.rq = rq1_1d;
end

function E = ghost_extrapolation(N, gp, p)
  Ng = N + 1 + 2 * p;
  Ni = N + 1 + 2 * gp;
  EI = (1:Ni)' + (p - gp);
  EJ = (1:Ni)';
  EV = ones(Ni,1);
  E = sparse(EI, EJ, EV, Ng, Ni);
  E * [-gp:N+gp]';
  r_src = (p-gp) + (-p:p+1);
  r_dst = (-p:p+1);
  I = lagrange_interpolation_matrix(r_src, r_dst);
  E(1:(2*p+2),1:(2*p+2)) = I;
  E(end + 1 - (1:(2*p+2)),end + 1 -(1:(2*p+2))) = I;
end
