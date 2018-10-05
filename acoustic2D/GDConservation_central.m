clear

alpha = 0;

% GD polynomial order is 2*p+1
for p = [3]
  % Number of GD ghost points
  for Ng = [0,p]
    % Watertight mesh
    [T_CM{p}, eng_CM{p}, m_pr_CM{p}, m_vx_CM{p}, m_vy_CM{p}] =  ...
      GDConservation_func(p, [20 15 15 20], 1, Ng, alpha);

    % discontinuous mesh
    [T_DM{p}, eng_DM{p}, m_pr_DM{p}, m_vx_DM{p}, m_vy_DM{p}] = ...
      GDConservation_func(p, [20 15 15 20], 0, Ng, alpha);

    % reference mesh
    [T_CF{p}, eng_CF{p}, m_pr_CF{p}, m_vx_CF{p}, m_vy_CF{p}] =  ...
      GDConservation_func(p, [20 20 20 20], 2, Ng, alpha);

    figure(1)
    semilogy(T_CM{p}, eng_CM{p}(2:end)/eng_CM{p}(1), '-', ...
             T_DM{p}, eng_DM{p}(2:end)/eng_DM{p}(1), '-', ...
             T_CF{p}, eng_CF{p}(2:end)/eng_CF{p}(1), '-')
    title('time versus energy')
    legend('watertight', 'discontinuous', 'reference')
    matlab2tikz(sprintf('GDbox_energy_n%d_Ng%d_central.tikz', 2*p+1, Ng))

    figure(2)
    semilogy(T_CM{p}, abs(m_pr_CM{p}(2:end) - m_pr_CM{p}(1)), '-', ...
             T_DM{p}, abs(m_pr_DM{p}(2:end) - m_pr_DM{p}(1)), '-', ...
             T_CF{p}, abs(m_pr_CF{p}(2:end) - m_pr_CF{p}(1)), '-')
    title('time versus conservation of pressure error')
    legend('watertight', 'discontinuous', 'reference')
    matlab2tikz(sprintf('GDbox_mass_pr_n%d_Ng%d_central.tikz', 2*p+1, Ng))

    figure(3)
    semilogy(T_CM{p}, abs(m_vx_CM{p}(2:end) - m_vx_CM{p}(1)), '-', ...
             T_DM{p}, abs(m_vx_DM{p}(2:end) - m_vx_DM{p}(1)), '-', ...
             T_CF{p}, abs(m_vx_CF{p}(2:end) - m_vx_CF{p}(1)), '-')
    title('time versus conservation of vx error')
    matlab2tikz(sprintf('GDbox_mass_vx_n%d_Ng%d_central.tikz', 2*p+1, Ng))

    figure(4)
    semilogy(T_CM{p}, abs(m_vy_CM{p}(2:end) - m_vy_CM{p}(1)), '-', ...
             T_DM{p}, abs(m_vy_DM{p}(2:end) - m_vy_DM{p}(1)), '-', ...
             T_CF{p}, abs(m_vy_CF{p}(2:end) - m_vy_CF{p}(1)), '-')
    title('time versus conservation of vy error')
    matlab2tikz(sprintf('GDbox_mass_vy_n%d_Ng%d_central.tikz', 2*p+1, Ng))
  end
end
