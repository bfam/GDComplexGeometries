function dt = max_Taylor_timestep(q, lam)
  lam(real(lam) > 0) = 0;

  dt = 1;
  v = max(abs(T(q, dt, lam)));

  dt = [dt, dt];
  if v <= 1
    while v < 1
      dt(2) = 2 * dt(2);
      v = max(abs(T(q, dt(2), lam)));
    end
  else
    while v > 1
      dt(1) = dt(1) / 2;
      v = max(abs(T(q, dt(1), lam)));
    end
  end

  tol = 1e-12;
  while diff(dt) > (1+max(dt)) * tol
    new_dt = (dt(1) + dt(2)) / 2;
    v = max(abs(T(q, new_dt, lam)));
    if v <= 1
      dt(1) = new_dt;
    else
      dt(2) = new_dt;
    end
  end

  plot_Taylor_stability_region(q)
  hold on
  plot(lam * dt(1), 'o')
  hold off

  dt = (dt(1) + dt(2)) / 2;


end

function e = T(q, dt, lam)

  e = ones(size(lam));
  de = e;

  for k = 1:q
    de = lam .* de * dt / k;
    e = e + de;
  end
end
