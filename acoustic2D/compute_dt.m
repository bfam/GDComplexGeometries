% estimate the time step for the problem defined by OP
function dt = compute_dt(OP)

Globals2D_gddg;

dt = inf;
for k = 1:length(OP.B)
  if ~isfield(OP.B{k},'isdg')
    dt = min(dt, 0.18 * (2 - OP.B{k}.Ng / OP.B{k}.p) * OP.B{k}.h);
  end
end

if isfield(OP, 'dg')
  Globals2D;
  rLGL = JacobiGQ(0,0,N); rmin = abs(rLGL(1)-rLGL(2));
  dtscale = dtscale2D;
  dt_dg =  min(dtscale)*rmin;
  % disp([dt dt_dg])
  dt = min(dt, dt_dg);
end
% fprintf('dt = %e\n', dt);
