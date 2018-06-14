% Compute the L2 error and energy of the operators (compute component integrals
% if compute_mass == true). See Acoustic2D for usage.
function [err, eng, err_gd, eng_gd, err_dg, eng_dg, m_pr, m_v1, m_v2] = ...
  L2Error2D(OP, time, v1, v2, pr, compute_mass)

Globals2D_gddg;

err_gd = 0;

%%%%%%%%%%%%%%%%%%%%%%%
% Galerkin Difference %
%%%%%%%%%%%%%%%%%%%%%%%

err_gd = 0;
eng_gd = 0;
m_pr = 0;
m_v1 = 0;
m_v2 = 0;

for k = 1:length(OP.B)

  if isfield(OP.B{k}, 'isdg')
    Globals2D;
    continue
  end

  B = OP.B{k};
  err = B.error;

  x1g = err.x1g;
  x2g = err.x2g;

  if isfield(B, 'E1')
    E = kron(B.E1, B.E2);
    Pv1 = err.P1F*(err.P2F*(E*v1(B.vmap(:))));
    Pv2 = err.P1F*(err.P2F*(E*v2(B.vmap(:))));
    Ppr = err.P1F*(err.P2F*(E*pr(B.vmap(:))));
  else
    Pv1 = err.P1F*(err.P2F*v1(B.vmap(:)));
    Pv2 = err.P1F*(err.P2F*v2(B.vmap(:)));
    Ppr = err.P1F*(err.P2F*pr(B.vmap(:)));
  end

  if isfield(B, 'isaffine')
    PJ = B.J;
  else
    PJ = err.J;
  end

  if isa(Exact2D, 'function_handle')
    [ev1, ev2, epr] = Exact2D(time, x1g, x2g);
  else
    ev1 = Pv1;
    ev2 = Pv2;
    epr = Ppr;
  end

  err_gd = err_gd + sum(err.w .* PJ .* (Pv1 - ev1).^2) / 2;
  err_gd = err_gd + sum(err.w .* PJ .* (Pv2 - ev2).^2) / 2;
  err_gd = err_gd + sum(err.w .* PJ .* (Ppr - epr).^2) / 2;

  eng_gd = eng_gd + sum(err.w .* PJ .* Pv1.^2) / 2;
  eng_gd = eng_gd + sum(err.w .* PJ .* Pv2.^2) / 2;
  eng_gd = eng_gd + sum(err.w .* PJ .* Ppr.^2) / 2;

  %{
  dv1 = B.P2F' * (B.P1F' * (B.w .* ((B.P1F * (B.P2F * v1(B.vmap(:)))) - ev1)));
  dv2 = B.P2F' * (B.P1F' * (B.w .* ((B.P1F * (B.P2F * v2(B.vmap(:)))) - ev2)));
  dpr = B.P2F' * (B.P1F' * (B.w .* ((B.P1F * (B.P2F * pr(B.vmap(:)))) - epr)));

  err_gd = err_gd + (dv1' * (B.MJI \ dv1));
  err_gd = err_gd + (dv2' * (B.MJI \ dv2));
  err_gd = err_gd + (dpr' * (B.MJI \ dpr));
  %}

  if compute_mass
      Pv1 = B.P2F' * (B.P1F' * (B.w .* (B.P1F * (B.P2F * v1(B.vmap(:))))));
      Pv2 = B.P2F' * (B.P1F' * (B.w .* (B.P1F * (B.P2F * v2(B.vmap(:))))));
      Ppr = B.P2F' * (B.P1F' * (B.w .* (B.P1F * (B.P2F * pr(B.vmap(:))))));
    m_pr = m_pr + sum(B.w .* (B.P1F*(B.P2F*((B.MJI \ Ppr)))));
    m_v1 = m_v1 + sum(B.w .* (B.P1F*(B.P2F*((B.MJI \ Pv1)))));
    m_v2 = m_v2 + sum(B.w .* (B.P1F*(B.P2F*((B.MJI \ Pv2)))));
  end

end

err = err_gd;
eng = eng_gd;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discontinuous Galerkin %
%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(OP, 'dg')

  if isfield(OP.dg, 'cubature_err')
    c = OP.dg.cubature_err;
  else
    c = OP.dg.cubature;
  end
  Pv1 = c.V*v1(vmapDG);
  Pv2 = c.V*v2(vmapDG);
  Ppr = c.V*pr(vmapDG);

  x1g = c.V*OP.x1(vmapDG);
  x2g = c.V*OP.x2(vmapDG);

  if isa(Exact2D, 'function_handle')
    [ev1, ev2, epr] = Exact2D(time, x1g, x2g);
  else
    ev1 = Pv1;
    ev2 = Pv2;
    epr = Ppr;
  end

  err_dg =          sum(sum(diag(c.w) * (c.J .* (Pv1 - ev1).^2))) / 2;
  err_dg = err_dg + sum(sum(diag(c.w) * (c.J .* (Pv2 - ev2).^2))) / 2;
  err_dg = err_dg + sum(sum(diag(c.w) * (c.J .* (Ppr - epr).^2))) / 2;

  eng_dg =          sum(sum(diag(c.w) * (c.J .* (Pv1).^2))) / 2;
  eng_dg = eng_dg + sum(sum(diag(c.w) * (c.J .* (Pv2).^2))) / 2;
  eng_dg = eng_dg + sum(sum(diag(c.w) * (c.J .* (Ppr).^2))) / 2;

  err = err_dg + err_gd;
  eng = eng_dg + eng_gd;

  if compute_mass
    m_pr = nan;
    m_v1 = nan;
    m_v2 = nan;
  end
else
  err_dg = nan;
  eng_dg = nan;
end

err = sqrt(err);
eng = sqrt(eng);
err_dg = sqrt(err_dg);
eng_dg = sqrt(eng_dg);
err_gd = sqrt(err_gd);
eng_gd = sqrt(eng_gd);