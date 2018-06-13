% Spectral derivative matrix for polynomials defined on grid points x
function D = spectral_derivative(x)
  n = length(x);
  D = zeros(n,n);

  w = lagrange_weights(x);

  for k = 1:n
    for j = 1:n
      if k == j
        for l = 1:n
          if l ~= k
            D(j, k) = D(j, k) + 1 / (x(k) - x(l));
          end
        end
      else
        D(j, k) = (w(k) / w(j)) / (x(j) - x(k));
      end
    end
  end
end
