function [P_fd2gq, w_gq, r_gd, r_gq, D_gq, P_fd2ge] = ...
          gd_quadrature(p, N_q, N, Np, N_plot_points)

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
  [r_q, w_q] = gauss_cofs(N_q+1);
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
  % r_gq = P_fd2gq*r_gd;
  D_gq = kron(speye(N), D_q);

  r_e = linspace(-1, 1, N_plot_points) / N;
  Np_e = length(r_e);

  % projection to the gauss points
  P_b2e = zeros(Np_e, Np_b);
  for k = 1:Np_e
    for j = 1:Np_b
      P_b2e(k, j) = w_b(j) / (r_e(k) - r_b(j));
      if ~isfinite(P_b2e(k,j))
        P_b2e(k,:) = 0;
        P_b2e(k,j) = 1;
        break
      end
    end
    d = sum(P_b2e(k, :));
    P_b2e(k, :) = P_b2e(k, :) / d;
  end

  % projection to the intermediate Gauss points
  P_fd2ge = zeros(Np_e * N, Np);
  r_ge = zeros(Np_e * N,1);
  i = 0;
  h = 2 / N;
  for k = 0:N-1
    P_fd2ge(k * Np_e + (1:Np_e), k + (1:Np_b)) = P_b2e;
    r_ge(k * Np_e + (1:Np_e)) = r_e + h / 2 + k * h - 1;
  end
  P_fd2ge = sparse(P_fd2ge);
end
