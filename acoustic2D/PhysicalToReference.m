% Determine the reference coordinates for a physical point
function [r, s] = PhysicalToReference(xout, yout, r, s, x, y,  ...
                                      xr, yr, xs, ys, invV, N)

tol = 100 * eps(1.0);

for iter = 0:1000
  I = Vandermonde2D(N, r, s) * invV;

  fx = xout - I * x;
  fy = yout - I * y;
  if iter == 1
  end

  if(abs(fx) < tol && abs(fy) < tol)
    return
  end

  Jxr = -I*xr;
  Jxs = -I*xs;
  Jyr = -I*yr;
  Jys = -I*ys;

  J = Jxr * Jys - Jxs * Jyr;

  dr = (Jys * fx - Jxs * fy) / J;
  ds = (Jxr * fy - Jyr * fx) / J;

  r = r - dr;
  s = s - ds;

end

error('PhysicalToReference did not converge')
