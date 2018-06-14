% Compute the matrix for dq/dt = A * q using the naive method of multiplying by
% vectors of all zeros except a single 1 (e.g., find one column of A at a time)
function A = compute(OP)

Np = numel(OP.x1);
A = eye(3*Np);
for k = 1:3*Np
  disp([k, 3*Np])
  v1 = A(1:Np,k);
  v2 = A(Np + (1:Np),k);
  pr = A(2*Np + (1:Np),k);
  [A(1:Np,k), A(Np + (1:Np),k), A(2*Np + (1:Np),k)] = ...
     AcousticRHS2D(OP, 0, v1, v2, pr);
end
