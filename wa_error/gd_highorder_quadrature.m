function [r, s, rq, sq, W, I, Dr, Ds] = ...
                gd_highorder_quadrature(Nr, Ns, p, Ng, Np_q)

  [Ir, wr, r, rq, Dr_gq] = gd_quadrature(p, Nr, Nr + 2 * p + 1, Np_q);
  [Is, ws, s, sq, Ds_gq] = gd_quadrature(p, Ns, Ns + 2 * p + 1, Np_q);
  Er = ghost_extrapolation(Nr, Ng, p);
  Es = ghost_extrapolation(Ns, Ng, p);
  r = r(p-Ng+1:end-(p-Ng));
  s = s(p-Ng+1:end-(p-Ng));
  Ir = Ir * Er;
  Is = Is * Es;

  [r, s] = meshgrid(r, s); r = r(:); s = s(:);
  [rq, sq] = meshgrid(rq, sq); rq = rq(:); sq = sq(:);

  I = kron(Ir, Is);
  W = kron(wr, ws);
  Dr = kron(Dr_gq, eye(size(Ds_gq)));
  Ds = kron(eye(size(Dr_gq)), Ds_gq);

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

function [P_fd2gq, w_gq, r_gd, r_gq, D_gq] = gd_quadrature(p, N, Np, Np_q)

  % Check the input arguments
  if(Np ~= N + 2 * p + 1)
    error('Np == N + 2 * p + 1 is required')
  end

  % Get the interpolation weights and points
  r_b = (2/N) * ((0:2*p+1)' - (p+0.5));
  w_b = zeros(size(r_b));
  for k = 0:2 * p + 1
    w_b(k + 1) = (-1)^k * nchoosek(2*p+1, k);
  end
  Np_b = length(r_b);

  % Get the integration points and weights
  [r_q, w_q] = gauss_cofs(Np_q);
  r_q = r_q / N;
  w_q = w_q / N;
  Np_q = length(r_q);
  D_q = spectral_derivative(r_q);

  % projection to the gauss points
  P_b2q = zeros(Np_q, Np_b);
  for k = 1:Np_q
    for j = 1:Np_b
      P_b2q(k, j) = w_b(j) / (r_q(k) - r_b(j));
    end
    d = sum(P_b2q(k, :));
    P_b2q(k, :) = P_b2q(k, :) / d;
  end

  % projection to the intermediate Gauss points
  P_fd2gq = zeros(Np_q * N, Np);
  r_gq = zeros(Np_q * N,1);
  i = 0;
  h = 2 / N;
  for k = 0:N-1
    P_fd2gq(k * Np_q + (1:Np_q), k + (1:Np_b)) = P_b2q;
    r_gq(k * Np_q + (1:Np_q)) = r_q + h / 2 + k * h - 1;
  end
  P_fd2gq = sparse(P_fd2gq);


  % stacked Gauss integration weights
  w_gq = kron(ones(N,1),w_q);

  r_gd = (-p:N+p)'*(2 / N) - 1;
  D_gq = kron(speye(N), D_q);
end
