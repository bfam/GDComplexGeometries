% Compute the projection errors
% Inputs:
%    beta :: Coordinate transform parameters
%    N1   :: N1+1 interior GD points are used in r1 direction
%    N2   :: N2+1 interior GD points are used in r2 direction
%    p    :: GD polynomial order is 2*p+1
%    Ng   :: Number of GD ghost points (0 <= Ng <= p)
%    Np_q :: Number of quadrature points for "exact" integration
%    k    :: mode number for the approximation
% Outputs:
%    errL2   :: L2 error for "exact" L2 projection
%    errWA   :: L2 error for "exact" WADG projection
%    errJnWA :: L2 error for inexact WADG projection
function [errL2, errWA, errJnWA] = wa_error_fun(beta, N1, N2, p, Ng, Np_q, k)

% Coordinate transform
x = @(r,s) r + beta .* cos(3 * pi * s / 2) .* cos(pi * r / 2);
y = @(r,s) s + beta .* sin(3 * pi * r / 2) .* cos(pi * s / 2);

% exact metric terms
xr = @(r,s) 1 - ( pi / 2) * beta * cos(3 * pi * s / 2) .* sin(pi * r / 2);
xs = @(r,s) -(3 * pi / 2) * beta * sin(3 * pi * s / 2) .* cos(pi * r / 2);
yr = @(r,s)  (3 * pi / 2) * beta * cos(3 * pi * r / 2) .* cos(pi * s / 2);
ys = @(r,s) 1 - ( pi / 2) * beta * sin(3 * pi * r / 2) .* sin(pi * s / 2);
J = @(r,s) xr(r,s).*ys(r,s) - xs(r,s).*yr(r,s);

% "exact" integration
[r, s, rq, sq, W, I, Dr, Ds] = gd_highorder_quadrature(N1, N2, p, Ng, Np_q);

% function to approximate
f = @(r,s) cos(k * x(r,s)) .* cos(k * y(r,s));

% "exact" integration of the function
MJf = I' * (W .* J(rq,sq) .* f(rq, sq));

% "exact" mass matrix
MJ  = I' * (diag(sparse(W .* J(rq, sq))) * I);
% "exact" WADG inverse mass matrix
MJI = I' * (diag(sparse(W ./ J(rq, sq))) * I);
% reference mass matrix
M   = I' * (diag(sparse(W             )) * I);

% "exact" L2 projection
gL2 = MJ \ MJf;

% "exact" WADG projection
gWA = M \ (MJI * (M \ MJf));

% L2 error for L2 and WADG
errL2 = sum(W .* J(rq, sq) .* (f(rq, sq) - I * gL2).^2);
errWA = sum(W .* J(rq, sq) .* (f(rq, sq) - I * gWA).^2);

% variational crimes integration
[r2, s2, rq2, sq2, W2, I2, Dr2, Ds2] = gd_highorder_quadrature(N1, N2, p, Ng, 2*p+2);
M2 = I2' * (diag(sparse(W2)) * I2);

% approximate metric terms, compute at quadrature nodes then project to GD space
x = I2 * (M2 \ (I2' * (W2 .* x(rq2, sq2))));
y = I2 * (M2 \ (I2' * (W2 .* y(rq2, sq2))));
xr = Dr2 * x;
xs = Ds2 * x;
yr = Dr2 * y;
ys = Ds2 * y;
Jn = xr.*ys - xs.*yr;
Jn = I2 * (M2 \ (I2' * (W2 .* Jn)));

% inexact WADG inverse mass matrix
MJnI = I2' * (diag(sparse(W2 ./ Jn))) * I2;

% inexact WADG projection
gJnWA = M \ (MJnI * (M \ MJf));

% inexact WADG error
errJnWA = sum(W .* J(rq, sq) .* (f(rq, sq) - I * gJnWA).^2);

% Check if the Jacobian determinant was negative
if min(J(rq,sq)) < 0 || min(Jn) < 0
  disp([min(J(rq,sq)), min(Jn)])
  if min(J(rq,sq)) < 0
    errJn = -errJn;
  end
  if min(Jn) < 0
    errJnWA = -errJnWA;
  end
end

