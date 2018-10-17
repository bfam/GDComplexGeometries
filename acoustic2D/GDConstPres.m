clear

% GD polynomial order is 2*p+1
for p = 3
  % Number of GD ghost points
  for Ng = [0,p]
    % Watertight mesh
    [T_CM{p}, err_CM{p}] = GDConstPres_func(p, [20 15 15 20], 1, Ng);

    % discontinuous mesh
    [T_DM{p}, err_DM{p}] = GDConstPres_func(p, [20 15 15 20], 0, Ng);

    % reference mesh
    [T_CF{p}, err_CF{p}] = GDConstPres_func(p, [20 20 20 20], 2, Ng);

    semilogy(T_CM{p}, err_CM{p}, '-', ...
             T_DM{p}, err_DM{p}, '-', ...
             T_CF{p}, err_CF{p}, '-')
    matlab2tikz(sprintf('constant_preserving_error_n%d_Ng%d.tikz', 2*p+1, Ng))
  end
end
