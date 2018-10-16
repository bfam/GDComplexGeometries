% Setup curved GD operator
% Inputs:
%   N1             :: Interior grid in the r1 direction is size N1+1
%   N2             :: Interior grid in the r2 direction is size N2+1
%   p              :: GD polynomial order is 2*p+1
%   quad_order     :: order of quadrature to use for subcells
%   grid           :: grid data structure should contain functions:
%                       grid.x1(r1, r2)
%                       grid.x2(r1, r2)
%                     which define the (x1,x2) coordinates for this block
%                     can also optionally contain an array of length 4
%                       grid.interp = [f1 f2 f3 f4]
%                     which
%                       -> fj == 0 will replace the L2 projected coordinate
%                          points along face j with the GD point interpolated
%                          values from functions grid.x1 and grid.x2
%                       -> If fj > 0 then values of coordinates along face j
%                          will be replaced with coordinate values defined by
%                          interpolating grid.x1 and grid.x2 along face j to a
%                          polynomial of degree 2*p+1; Chebyshev points of the
%                          2nd kind are used
%                       -> If fj < 0 then values along face j remain unchanged
%                          (e.g., just the values that arose from the L2
%                          projection of the volume coordinate transform into
%                          the GD space)
%   Ng             :: (optional) Number of free ghost points [default Ng = p]
%   err_quad_order :: (optional) order of quadrature to use in the error
%                     calculation [default err_quad_order = quad_order]
%   compute_MJI_R  :: (optional) boolean for specifying whether Cholesky
%                     factorization of the weight-adjusted mass matrix should be
%                     computed (if computed this will be used for the energy
%                     calculation) [default compute_MJI_R = false]
% Outputs:
%   B :: GD operator struct
% members of the struct:
%   N1, N2                  :: interior grid dimensions
%   p                       :: polynomial order is 2*p+1
%   Ng                      :: number of ghost points
%   Np1, Np2                :: total number of points in each dimension
%   quad_order              :: subcell quadrature order
%   P1_ge, P2_ge            :: projections to equally spaced grid
%   E1, E2                  :: non-free ghost extrapolations
%   P1, P2                  :: interpolation to the quadrature grid
%   P2F, P1F                :: 2D quadrature interpolation [apply: P1F*(P2F*q)]
%   P2FT, P1FT              :: Transposes of P2F and P1F
%   w1, w2                  :: 1D quadrature grids
%   w                       :: 2D quadrature grid
%   r1_1d, r2_1d            :: 1D GD grids
%   rq1_1d, rq2_1d          :: 1D quadrature grids
%   M1g, M2g                :: 1D GD mass matrices
%   Dq1, Dq2                :: 1D quadrature grid differentiation matrices
%   rq1, rq2                :: 2D quadrature grid
%   r1, r2                  :: 2D GD grid
%   x1, x2                  :: 2D GD coordinate points
%   M                       :: 2D GD mass matrix
%   R1, R2                  :: Dimension by dimension Cholesky factors
%   D1, D2                  :: 2D differentiation matrices
%   x1_1, x2_1, x1_2, x2_2  :: metric derivatives (interpolated to quadrature
%                              grid)
%   J                       :: Jacobian determinant (interpolated to quadrature
%                              grid)
%   r1_1, r2_1, r1_2, r2_2  :: metric derivatives (interpolated to quadrature
%                              grid)
%   h                       :: measure of minimum grid spacing in physical
%                              coordinates
%   Sx, Sy                  :: functions to apply stiffness matrices
%   SxT, SyT                :: functions to apply stiffness matrices transpose
%   massInv                 :: function to apply WA mass matrix inverse
%   MJI                     :: function to apply 1/J weighted mass matrix
%   R1_1d, R2_1d            :: 1D mass matrix Cholesky factors
%   massJI                  :: function to multiply 1/J at quadrature nodes
%   massInv_noJ             :: function to apply reference mass inverse
%   f                       :: face map struct
%     -> vmap               :: volume points on this face
%     -> P                  :: interpolation to face quadrature
%     -> w                  :: face quadrature weights
%     -> r1, r2             :: face GD points
%     -> g                  :: face cell edges
%     -> rq                 :: face quadrature points
%     -> corners            :: corners of the face
%     -> x1, x2             :: physical coordinates of face GD points
%     -> sJ                 :: surface Jacobian
%     -> n1, n2             :: surface unit normal
%   error                   :: structure for error calculation
%     -> P2F, P1F           :: error quadrature interpolations
%     -> x1g, x2g           :: physicaly coordinates of error quadrature points
%     -> J                  :: Jacobian determinant at error quadrature points
%     -> w                  :: error quadrature weights
function [B] = gd_setup_curved_alias(N1, N2, p, quad_order, grid, Ng, ...
                                     err_quad_order, compute_MJI_R)

  if nargin < 6
    Ng = [];
  end
  if isempty(Ng)
    Ng = p;
  end

  if nargin < 7
    err_quad_order = [];
  end
  if isempty(err_quad_order)
    err_quad_order = quad_order;
  end

  if nargin < 8
    compute_MJI_R = [];
  end
  if isempty(compute_MJI_R)
    compute_MJI_R = false;
  end


  N_plot_points = 10;

  B.N1  = N1;
  B.N2  = N2;
  B.p   = p;
  B.Ng  = p;

  % number of points with ghost
  Ng1 = B.N1 + 2 * B.p + 1;
  Ng2 = B.N2 + 2 * B.p + 1;
  Np1 = Ng1 - 2 * (p - Ng);
  Np2 = Ng2 - 2 * (p - Ng);
  B.Np1 = Np1;
  B.Np2 = Np2;

  % compute the gd quadrature rule (all operators assume free ghost points)
  B.quad_order = quad_order;
  % P1    :: interpolation from GD to quadrature grid
  % w1    :: quadrature weights
  % r1    :: GD grid (with ghost points)
  % rq1   :: quadrature grid
  % Dq1   :: differentiation matrix defined for quadrature grid values
  % P1_ge :: interpolation to equally spaced grid (useful for plotting)
  [P1, w1, r1, rq1, Dq1, B.P1_ge] = ...
    gd_quadrature(p, B.quad_order, N1, Ng1, N_plot_points);
  [P2, w2, r2, rq2, Dq2, B.P2_ge] = ...
    gd_quadrature(p, B.quad_order, N2, Ng2, N_plot_points);

  % Build operators to extrapolate the ghost values
  E1 = ghost_extrapolation(N1, Ng, p);
  E2 = ghost_extrapolation(N2, Ng, p);
  B.E1 = E1;
  B.E2 = E2;

  % Chain the ghost operators with the free operators to get the real operators
  B.P1_ge = B.P1_ge * E1;
  B.P2_ge = B.P2_ge * E2;
  P1 = P1 * E1;
  P2 = P2 * E2;
  B.P1 = P1;
  B.P2 = P2;

  % GD grid of degrees of freedom (remove non free ghost points)
  r1 = r1((p-Ng+1):(N1+1+2*p-p+Ng));
  r2 = r2((p-Ng+1):(N2+1+2*p-p+Ng));

  % 2-D interpolation operators
  % (need to be applied in order B.P1F * (B.P2F * q))
  B.P2F = kron(speye(size(B.P1,2)), P2);
  B.P1F = kron(P1, speye(size(B.P2,1)));
  B.P2FT = B.P2F';
  B.P1FT = B.P1F';

  % 2-D quadrature weights
  B.w1 = w1;
  B.w2 = w2;
  B.w  = kron(w1, w2);

  B.r1_1d = r1;
  B.r2_1d = r2;

  B.rq1_1d = rq1;
  B.rq2_1d = rq2;

  % 1-D mass matrices (for GD based functions)
  B.M1g = P1' * diag(sparse(w1)) * P1;
  B.M2g = P2' * diag(sparse(w2)) * P2;

  R1 = chol(B.M1g);
  R2 = chol(B.M2g);

  B.Dq1 = Dq1;
  B.Dq2 = Dq2;

  % Set up reference grid
  [rq1,rq2] = meshgrid(rq1, rq2);
  B.rq1 = rq1(:);
  B.rq2 = rq2(:);
  xq1 = grid.x1(rq1, rq2);
  xq2 = grid.x2(rq1, rq2);

  % L2 projection of (x,y) to GD basis
  x1 = R2 \ (R2' \ (P2' * (((((w2 * w1') .* xq1) * P1) / R1) / R1')));
  x2 = R2 \ (R2' \ (P2' * (((((w2 * w1') .* xq2) * P1) / R1) / R1')));

  % 2-D GD reference mesh
  r1_1d = r1;
  r2_1d = r2;
  [r1,r2] = meshgrid(r1, r2);
  B.r1 = r1(:);
  B.r2 = r2(:);

  % Project out the face modes is necessary
  if isfield(grid, 'interp')
    q = 2*p+1;
    if grid.interp(1) > 0
      x1(:,Ng+1) = polynomial_interpolate_fn(r2_1d, @(s) grid.x1(-1, s), q);
      x2(:,Ng+1) = polynomial_interpolate_fn(r2_1d, @(s) grid.x2(-1, s), q);
    elseif grid.interp(1) == 0
      x1(:,Ng+1) = grid.x1(-1, r2_1d);
      x2(:,Ng+1) = grid.x2(-1, r2_1d);
    end
    if grid.interp(2) > 0
      x1(:,Ng+N1+1) = polynomial_interpolate_fn(r2_1d, @(s) grid.x1(1, s), q);
      x2(:,Ng+N1+1) = polynomial_interpolate_fn(r2_1d, @(s) grid.x2(1, s), q);
    elseif grid.interp(2) == 0
      x1(:,Ng+N1+1) = grid.x1(1, r2_1d);
      x2(:,Ng+N1+1) = grid.x2(1, r2_1d);
    end
    if grid.interp(3) > 0
      x1(Ng+1,:) = polynomial_interpolate_fn(r1_1d, @(r) grid.x1(r, -1), q);
      x2(Ng+1,:) = polynomial_interpolate_fn(r1_1d, @(r) grid.x2(r, -1), q);
    elseif grid.interp(3) == 0
      x1(Ng+1, :) = grid.x1(r1_1d, -1);
      x2(Ng+1, :) = grid.x2(r1_1d, -1);
    end
    if grid.interp(4) > 0
      x1(Ng+N2+1,:) = polynomial_interpolate_fn(r1_1d, @(r) grid.x1(r, 1), q);
      x2(Ng+N2+1,:) = polynomial_interpolate_fn(r1_1d, @(r) grid.x2(r, 1), q);
    elseif grid.interp(4) == 0
      x1(Ng+N2+1, :) = grid.x1(r1_1d, 1);
      x2(Ng+N2+1, :) = grid.x2(r1_1d, 1);
    end
  end

  B.x1 = x1(:);
  B.x2 = x2(:);

  % 2-D mass matrix
  B.M   = kron(       B.M1g,        B.M2g);

  % TODO: not needed?
  %% B.M2  = kron(speye(Np1),       B.M2g);
  %% B.M1  = kron(       B.M1g,  speye(Np2));


  % dimension by dimension Cholesky factors
  B.R1  = kron( chol(B.M1g),  speye(Np2));
  B.R2  = kron(speye(Np1),  chol(B.M2g));
  % B.R   = kron(chol(B.M1g), chol(B.M2g));

  % 2-D differentiation matrices
  B.D1  = kron(B.Dq1, speye(length(w2)));
  B.D2  = kron(speye(length(w1)), B.Dq2);

  %
  % Setup the metric terms
  %

  % Compute metric derivatives and project back
  Dr = kron(R1 \ (R1' \ (P1' * diag(sparse(w1)) * B.Dq1 * P1)),speye(Np2));
  Ds = kron(speye(Np1),R2 \ (R2' \ (P2' * diag(sparse(w2)) * B.Dq2 * P2)));
  B.x1_1 = Dr * B.x1;
  B.x2_1 = Dr * B.x2;
  B.x1_2 = Ds * B.x1;
  B.x2_2 = Ds * B.x2;

  xg1_1 = B.P1F * (B.P2F * B.x1_1);
  xg2_1 = B.P1F * (B.P2F * B.x2_1);
  xg1_2 = B.P1F * (B.P2F * B.x1_2);
  xg2_2 = B.P1F * (B.P2F * B.x2_2);

  B.J = xg1_1 .* xg2_2 - xg1_2 .* xg2_1;

  B.r1_1 =  xg2_2 ./ B.J;
  B.r2_1 = -xg2_1 ./ B.J;
  B.r1_2 = -xg1_2 ./ B.J;
  B.r2_2 =  xg1_1 ./ B.J;

  h1 = 2 * min(sqrt(B.x1_1.^2 + B.x2_1.^2))/N1;
  h2 = 2 * min(sqrt(B.x1_2.^2 + B.x2_2.^2))/N2;
  B.h = min(h1,h2);

  S1 = @(m1, m2) B.P2FT * (B.P1FT * (diag(sparse(B.w .* m1 .* m2)) * B.D1) * B.P1F) * B.P2F;
  S2 = @(m1, m2) B.P2FT * (B.P1FT * (diag(sparse(B.w .* m1 .* m2)) * B.D2) * B.P1F) * B.P2F;
  Ax = S1(B.J, B.r1_1) + S2(B.J, B.r2_1);
  Ay = S1(B.J, B.r1_2) + S2(B.J, B.r2_2);

  B.Sx = @(v) Ax * v(:);
  B.Sy = @(v) Ay * v(:);
  B.SxT = @(v) Ax' * v(:);
  B.SyT = @(v) Ay' * v(:);

  if min(B.J(:)) < 0
     error('negative Jacobian')
  end

  % to simplify our lives in Matlab, we build the matrix
  MJI = B.P2FT * (B.P1FT * diag(sparse(B.w ./ B.J)) * B.P1F) * B.P2F;
  B.massInv = @(v) (B.R1 \ (B.R1' \ (B.R2 \ (B.R2' \ (MJI * (B.R2 \ (B.R2' \ (B.R1 \ (B.R1' \ v(:))))))))));
  B.MJI = MJI;
  if compute_MJI_R
    B.MJI_R = chol(MJI);
  end

  vec_wJI = reshape(B.w ./ B.J, numel(B.w2), numel(B.w1));

  B.R1_1d  = chol(B.M1g);
  B.R2_1d  = chol(B.M2g);
  B.massJI = @(v) reshape(MJI * v(:), size(v));

  B.massInv_noJ = @(v) (B.R2_1d \ (B.R2_1d' \ ((v / B.R1_1d) / B.R1_1d')));
  B.massInv = @(v) B.massInv_noJ(B.massJI(B.massInv_noJ(v)));

  %
  % Face stuff
  %

  % Boundary maps
  B.f{1}.vmap  =  Ng      *Np2 +       (1:Np2)';
  B.f{2}.vmap  = (Ng+B.N1)*Np2 +       (1:Np2)';
  B.f{3}.vmap = Ng+1           + Np2*(0:Np1-1)';
  B.f{4}.vmap = Ng+B.N2+1      + Np2*(0:Np1-1)';

  B.f{1}.P = B.P2; B.f{1}.w = B.w2;
  B.f{2}.P = B.P2; B.f{2}.w = B.w2;
  B.f{3}.P = B.P1; B.f{3}.w = B.w1;
  B.f{4}.P = B.P1; B.f{4}.w = B.w1;

  %
  % Boundary Metric Terms
  %
  B.f{1}.r1 = -ones(size(B.r2_1d));                  % gd space r1
  B.f{1}.r2 = B.r2_1d;                               % gd space r2
  B.f{1}.g  = linspace(-1,1,B.N2+1)';                % gd cells edges
  B.f{1}.rq = B.rq2_1d;                              % quadrature points
  % B.f{1}.x1 = grid.x1(B.f{1}.r1, B.f{1}.r2);         % gd space x1
  % B.f{1}.x2 = grid.x2(B.f{1}.r1, B.f{1}.r2);         % gd space x2
  B.f{1}.corners = [grid.x1(-1,-1), grid.x2(-1,-1);  % bottom corner of the edge
                    grid.x1(-1, 1), grid.x2(-1, 1)]; % bottom corner of the edge
  B.f{1}.x1 = x1(B.f{1}.vmap);
  B.f{1}.x2 = x2(B.f{1}.vmap);
  x1_2 = B.f{1}.P * B.x1_2(B.f{1}.vmap);             % dx1/dr in quadrature space
  x2_2 = B.f{1}.P * B.x2_2(B.f{1}.vmap);             % dx2/dr in quadrature space
  B.f{1}.sJ =  sqrt(x1_2.^2 + x2_2.^2);              % surface Jacobian at quadrature nodes
  B.f{1}.n1 = -x2_2 ./ B.f{1}.sJ;                    % n1 at quadrature nodes
  B.f{1}.n2 =  x1_2 ./ B.f{1}.sJ;                    % n2 at quadrature nodes

  B.f{2}.r1 = ones(size(B.r2_1d));
  B.f{2}.r2 = B.r2_1d;
  B.f{2}.g  = linspace(-1,1,B.N2+1)';
  B.f{2}.rq = B.rq2_1d;
  % B.f{2}.x1 = grid.x1(B.f{2}.r1, B.f{2}.r2);
  % B.f{2}.x2 = grid.x2(B.f{2}.r1, B.f{2}.r2);
  B.f{2}.corners = [grid.x1(1,-1), grid.x2(1,-1);
                    grid.x1(1, 1), grid.x2(1, 1)];
  B.f{2}.x1 = x1(B.f{2}.vmap);
  B.f{2}.x2 = x2(B.f{2}.vmap);
  x1_2 = B.f{2}.P * B.x1_2(B.f{2}.vmap);
  x2_2 = B.f{2}.P * B.x2_2(B.f{2}.vmap);
  B.f{2}.sJ =  sqrt(x1_2.^2 + x2_2.^2);
  B.f{2}.n1 =  x2_2 ./ B.f{2}.sJ;
  B.f{2}.n2 = -x1_2 ./ B.f{2}.sJ;

  B.f{3}.r1 = B.r1_1d;
  B.f{3}.r2 = -ones(size(B.r1_1d));
  B.f{3}.g  = linspace(-1,1,B.N1+1)';
  B.f{3}.rq = B.rq1_1d;
  % B.f{3}.x1 = grid.x1(B.f{3}.r1, B.f{3}.r2);
  % B.f{3}.x2 = grid.x2(B.f{3}.r1, B.f{3}.r2);
  B.f{3}.corners = [grid.x1(-1,-1), grid.x2(-1,-1);
                    grid.x1( 1,-1), grid.x2( 1,-1)];
  B.f{3}.x1 = x1(B.f{3}.vmap);
  B.f{3}.x2 = x2(B.f{3}.vmap);
  x1_1 = B.f{3}.P * B.x1_1(B.f{3}.vmap);
  x2_1 = B.f{3}.P * B.x2_1(B.f{3}.vmap);
  B.f{3}.sJ =  sqrt(x1_1.^2 + x2_1.^2);
  B.f{3}.n1 =  x2_1 ./ B.f{3}.sJ;
  B.f{3}.n2 = -x1_1 ./ B.f{3}.sJ;

  B.f{4}.r1 = B.r1_1d;
  B.f{4}.r2 = ones(size(B.r1_1d));
  B.f{4}.g  = linspace(-1,1,B.N1+1)';
  B.f{4}.rq = B.rq1_1d;
  % B.f{4}.x1 = grid.x1(B.f{4}.r1, B.f{4}.r2);
  % B.f{4}.x2 = grid.x2(B.f{4}.r1, B.f{4}.r2);
  B.f{4}.corners = [grid.x1(-1, 1), grid.x2(-1, 1);
                    grid.x1( 1, 1), grid.x2( 1, 1)];
  B.f{4}.x1 = x1(B.f{4}.vmap);
  B.f{4}.x2 = x2(B.f{4}.vmap);
  x1_1 = B.f{4}.P * B.x1_1(B.f{4}.vmap);
  x2_1 = B.f{4}.P * B.x2_1(B.f{4}.vmap);
  B.f{4}.sJ =  sqrt(x1_1.^2 + x2_1.^2);
  B.f{4}.n1 = -x2_1 ./ B.f{4}.sJ;
  B.f{4}.n2 =  x1_1 ./ B.f{4}.sJ;

  % Error calculation
  [P1, w1, r1, rq1, Dq1] = ...
    gd_quadrature(p, err_quad_order, N1, Ng1, N_plot_points);
  [P2, w2, r2, rq2, Dq2] = ...
    gd_quadrature(p, err_quad_order, N2, Ng2, N_plot_points);

  [rq1,rq2] = meshgrid(rq1, rq2);

  x1g = grid.x1(rq1, rq2);
  x2g = grid.x2(rq1, rq2);

  E.P2F = kron(speye(size(P1,2)), P2);
  E.P1F = kron(P1, speye(size(P2,1)));

  x1_1 = x1g * Dq1';
  x2_1 = x2g * Dq1';
  x1_2 = Dq2 * x1g;
  x2_2 = Dq2 * x2g;

  J = x1_1 .* x2_2 - x1_2 .* x2_1;

  E.x1g = x1g(:);
  E.x2g = x2g(:);
  E.J = J(:);
  E.w  = kron(w1, w2);

  B.error = E;
end

% interpolates the function x with an enpoint preserving polynomial of degree Nc
% Evaluates the resulting polynomial at points r_gd
function [xf] = polynomial_interpolate_fn(r_gd, x, Nc)
  % Chebyshev points of the 2nd kind
  rc_c = sin(pi * linspace(-1/2, 1/2, Nc+1)');
  Ic2gd = lagrange_interpolation_matrix(rc_c, r_gd);
  xf = Ic2gd * x(rc_c);
end

% Build the matrix which extrapolates to ghost points using N+1 GD interior
% points polynomial order 2*p+1 and gp free ghost points (if gp = p then E is
% the identity matrix)
function E = ghost_extrapolation(N, gp, p)
  Ng = N + 1 + 2 * p;   Ni = N + 1 + 2 * gp;
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
