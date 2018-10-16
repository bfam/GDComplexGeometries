% Run the Acoustic2D problem
% Inputs:
%   OP           :: Cell array of GD operators (with possible DG dummy cell
%   v1           :: initial condition for v1
%   v2           :: initial condition for v2
%   pr           :: initial condition for pr
%   dt           :: time step
%   Ntsteps      :: Number of time steps
%   taylor_order :: order of the Taylor time
%                   if empty use RK
%                   > 0 use next highest Taylor order that includes imaginary
%                       axis
%                   < 0 use next highest Taylor order -taylor_order
%   compute_mass :: (optional) compute the mass time history (default: false)
% Ouputs:
%   v1     :: v1 at final time
%   v2     :: v2 at final time
%   pr     :: pr at final time
%   T      :: time for err, eng, etc.
%   err    :: error at times T
%   eng    :: energy at time T
%   err_gd :: error for gd cell at times T
%   err_dg :: error for dg cell at times T
%   m_pr   :: integral of pr at times T (only filled if compute_mass == true)
%   m_v1   :: integral of v1 at times T (only filled if compute_mass == true)
%   m_v2   :: integral of v2 at times T (only filled if compute_mass == true)

function [v1,v2,pr,T,err,eng,err_gd,err_dg,m_pr,m_v1,m_v2] = ...
  Acoustic2D(OP, v1, v2, pr, dt, Ntsteps, taylor_order, compute_mass)

Globals2D_gddg;

if nargin < 8
  compute_mass = [];
end
if isempty(compute_mass)
  compute_mass = false;
end

% Runge-Kutta residual storage
resv1 = zeros(size(v1)); resv2 = zeros(size(v2)); respr = zeros(size(pr));

time = 0; tstep = 1; d = 0;

[~, eng0, ~, eng_gd0, ~, eng_dg0, m_pr0, m_v10, m_v20] = ...
  L2Error2D(OP, 0, v1, v2, pr, compute_mass);
fnum = 0;
m_pr = [];
m_v1 = [];
m_v2 = [];
tic
for tstep = 1:Ntsteps
  if isempty(taylor_order)
    for INTRK = 1:5   % inner multi-stage Runge-Kutta loop
      localtime = time + rk4c(INTRK)*dt;
      [dv1, dv2, dpr] = AcousticRHS2D(OP, localtime, v1, v2, pr);

      % initiate, increment Runge-Kutta residuals and update fields
      resv1 = rk4a(INTRK)*resv1 + dt*dv1;   v1 = v1+rk4b(INTRK)*resv1;
      resv2 = rk4a(INTRK)*resv2 + dt*dv2;   v2 = v2+rk4b(INTRK)*resv2;
      respr = rk4a(INTRK)*respr + dt*dpr;   pr = pr+rk4b(INTRK)*respr;
    end
  else
    dv1 = v1;
    dv2 = v2;
    dpr = pr;
    % 3, 4, 7, 8, 11, 12
    if taylor_order < 0
      taylor_order = -taylor_order;
    else
      if taylor_order == 1 || taylor_order == 2
        taylor_order = 3;
      elseif taylor_order == 5 || taylor_order == 6
        taylor_order = 7;
      elseif taylor_order == 9 || taylor_order == 10
        taylor_order = 11;
      elseif taylor_order == 13 || taylor_order == 14
        taylor_order = 15;
      end
    end
    for k = 1:taylor_order
      [dv1, dv2, dpr] = AcousticRHS2D(OP, inf, dv1, dv2, dpr);
      dv1 = dt * dv1 / k;
      dv2 = dt * dv2 / k;
      dpr = dt * dpr / k;
      v1 = v1 + dv1;
      v2 = v2 + dv2;
      pr = pr + dpr;
    end
  end
  step_time = toc / tstep;

  time = tstep * dt;

  if (mod(tstep, 10) == 0 || tstep == Ntsteps)
    fnum = fnum + 1;
    fprintf('time = %e (step = %4d)\n', time, tstep);

    d = d+1;
    [err(d), eng(d), err_gd(d), eng_gd(d), err_dg(d), eng_dg(d), m_pr(d), m_v1(d), m_v2(d)] = ...
      L2Error2D(OP, dt * tstep, v1, v2, pr, compute_mass);
    T(d) = dt * tstep;

    if isfield(OP, 'dg')
      fprintf('energy0 (gd, dg) = %e (%e, %e)\n', ...
        eng0, eng_gd0, eng_dg0)
      fprintf('energy  (gd, dg) = %e (%e, %e)\n', ...
        eng(end), eng_gd(end), eng_dg(end))
      fprintf('error   (gd, dg) = %e (%e, %e)\n', ...
        err(end), err_gd(end), err_dg(end))
    else
      fprintf('>> energy0 = %e\n', eng0);
      fprintf('>> energy  = %e\n', eng(end));
      fprintf('>> en/en0  = %e\n', eng(end) / eng0);
      fprintf('>> error   = %e\n', err(end));
      if compute_mass
        fprintf('>> dm_pr   = %+e\n', m_pr(end) - m_pr0);
        fprintf('>> dm_v1   = %+e\n', m_v1(end) - m_v10);
        fprintf('>> dm_v2   = %+e\n', m_v2(end) - m_v20);
      end
    end
    fprintf('average time step timing = %.2e s\n', step_time);
    fprintf('estimated time remaining = %.2e s (%.2e min) (%.2e hr)\n', ...
            step_time * (Ntsteps - tstep), ...
            step_time * (Ntsteps - tstep) / 60, ...
            step_time * (Ntsteps - tstep) / 60 / 60)
    if 2*eng0 < eng(end)
      return
    end
  end
end

fprintf('\n\n\n')
fprintf('FINAL time = %e (step = %4d)\n', time, tstep);
if isfield(OP, 'dg')
  fprintf('FINAL energy0 (gd, dg) = %e (%e, %e)\n', ...
    eng0, eng_gd0, eng_dg0)
  fprintf('FINAL energy  (gd, dg) = %e (%e, %e)\n', ...
    eng(end), eng_gd(end), eng_dg(end))
  fprintf('FINAL error   (gd, dg) = %e (%e, %e)\n', ...
    err(end), err_gd(end), err_dg(end))
else
  fprintf('FINAL energy0 = %e\n', eng0);
  fprintf('FINAL energy  = %e\n', eng(end));
  fprintf('FINAL error   = %e\n', err(end));
  if compute_mass
    fprintf('FINAL error   = %e\n', err(end));
    fprintf('FINAL dm_pr   = %+e\n', m_pr(end) - m_pr0);
    fprintf('FINAL dm_v1   = %+e\n', m_v1(end) - m_v10);
    fprintf('FINAL dm_v2   = %+e\n', m_v2(end) - m_v20);
    eng  = [eng0, eng];
    m_pr = [m_pr0, m_pr];
    m_v1 = [m_v10, m_v1];
    m_v2 = [m_v20, m_v2];
  end
end
