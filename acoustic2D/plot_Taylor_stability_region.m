% Plot stability boundary of degree q Taylor Series
% Inputs:
%   q :: Taylor time stepping order
function plot_Taylor_stability_region(q)

phi = exp(1i * linspace(0, 2 * q * pi, q * 1000 + 1)');

lam = zeros(size(phi));

max_iterations = 100;
tol=10^(-12);

for k=2:length(phi)

  % Initial guess
  lnow = lam(k - 1) + (phi(k) - phi(k - 1)) / (q * T(lam(k - 1), q - 1));

  for it = 1:max_iterations
    dlam = -(T(lnow,q) - phi(k)) / T(lnow, q - 1);
    lnow = lnow + dlam;
    if dlam < tol
      break
    end
  end

  if (dlam > tol)
    error('failure k = %d', k)
  else
    lam(k)=lnow;
  end
end

plot(lam, '-')
title(sprintf('Taylor stability domain: order = %d',q))

function eT = T(l, q)

eT = 1;
df = l;

for j=1:q
  eT = eT + df;
  df = df * l / (j + 1);
end
