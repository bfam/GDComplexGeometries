% RHS evaluation of the acoustic wave equation
function [rhsv1, rhsv2, rhspr] = Acoustic2D(OP, time, v1, v2, pr)

Globals2D_gddg;

rhsv1 = zeros(size(v1)); rhsv2 = zeros(size(v2)); rhspr = zeros(size(pr));

alpha = 1;

for k = 1:length(OP.B)
  B = OP.B{k};

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Galerkin Difference Volume %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if ~isfield(B,'isdg')
    rhsv1(B.vmap) = B.SxT(pr(B.vmap));
    rhsv2(B.vmap) = B.SyT(pr(B.vmap));
    rhspr(B.vmap) = -(B.Sx(v1(B.vmap)) + B.Sy(v2(B.vmap)));
  end


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Galerkin Difference Faces and DG Mortar %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % faces
  km = k;
  for fm = 1:length(OP.B{km}.toB)
    % Get the plus side info
    kp = OP.B{km}.toB(fm);
    fp = OP.B{km}.toF(fm);

    % minus size
    B_vmapM = B.vmap(B.f{fm}.vmap);
    prM = B.f{fm}.P * pr(B_vmapM);
    v1M = B.f{fm}.P * v1(B_vmapM);
    v2M = B.f{fm}.P * v2(B_vmapM);
    vnM = B.f{fm}.n1 .* v1M + B.f{fm}.n2 .* v2M;

    % plus side
    if kp == 0
      % bc
      % prP = -prM; vnP = -vnM; % Zero pressure
      prP = prM; vnP = vnM; % Zero velocity
    else
      % plus from another block
      afp = abs(fp);
      B_vmapP = OP.B{kp}.vmap(OP.B{kp}.f{afp}.vmap);
      prP = OP.B{kp}.f{afp}.P * pr(B_vmapP);
      v1P = OP.B{kp}.f{afp}.P * v1(B_vmapP);
      v2P = OP.B{kp}.f{afp}.P * v2(B_vmapP);
      vnP = OP.B{kp}.f{afp}.n1 .* v1P + OP.B{kp}.f{afp}.n2 .* v2P;
    end

    prS = (prM + prP) / 2 + alpha * (vnM + vnP) / 2;
    vnS = (vnM - vnP) / 2 + alpha * (prM - prP) / 2;

    rhspr(B_vmapM) = rhspr(B_vmapM) + ...
      B.f{fm}.P' * (B.f{fm}.sJ .* B.f{fm}.w .* (vnM - vnS));
    rhsv1(B_vmapM) = rhsv1(B_vmapM) - ...
      B.f{fm}.P' * (B.f{fm}.sJ .* B.f{fm}.w .* B.f{fm}.n1 .* prS);
    rhsv2(B_vmapM) = rhsv2(B_vmapM) - ...
      B.f{fm}.P' * (B.f{fm}.sJ .* B.f{fm}.w .* B.f{fm}.n2 .* prS);

  end

  % multiply through inverse volume mass matrix and Jacobian
  if ~isfield(B, 'isdg')
    rhspr(B.vmap) = B.massInv(rhspr(B.vmap));
    rhsv1(B.vmap) = B.massInv(rhsv1(B.vmap));
    rhsv2(B.vmap) = B.massInv(rhsv2(B.vmap));
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discontinuous Galerkin Interior and Boundary %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(OP,'dg')
  Globals2D;

  %%%%%%%%%%%%%%%%%%
  % Surface Update %
  %%%%%%%%%%%%%%%%%%
  % g.W is Gauss weights and surface Jacobian
  g = OP.dg.gauss;

  % storage for field differences at faces
  gpr = g.interp * pr(vmapDG);
  gvn = g.nx .* (g.interp * v1(vmapDG)) + g.ny .* (g.interp * v2(vmapDG));

  gprM = gpr(g.mapM); gprP = gpr(g.mapP);
  gvnM = gvn(g.mapM); gvnP = gvn(g.mapP);

  % Fix the boundary

  % gprP(g.mapB) = -gprM(g.mapB); gvnP(g.mapB) = -gvnM(g.mapB); % Zero pressure
  gprP(g.mapB) =  gprM(g.mapB); gvnP(g.mapB) =  gvnM(g.mapB); % Zero velocity

  % Exact BC
  % [v1B, v2B, prB] = Exact2D(time, g.x, g.y);
  % vnB = g.nx .* v1B + g.ny .* v2B;
  % gprP(g.mapB) = prB(g.mapB); gvnP(g.mapB) = -vnB(g.mapB);

  % Zero out the mortar
  gprM(g.mapMor) = 0; gvnM(g.mapMor) = 0;
  gprP(g.mapMor) = 0; gvnP(g.mapMor) = 0;
  % [v1B, v2B, prB] = Exact2D(time, g.x, g.y);
  % gprP(g.mapMor) = prB(g.mapMor); gvnP(g.mapMor) = -vnB(g.mapMor);

  % Compute the numerical flux state
  gprS = (gprM + gprP) / 2 + alpha * (gvnM + gvnP) / 2;
  gvnS = (gvnM - gvnP) / 2 + alpha * (gprM - gprP) / 2;

  rhsv1(vmapDG) = rhsv1(vmapDG) - g.interp' * (g.W .* g.nx .* gprS);
  rhsv2(vmapDG) = rhsv2(vmapDG) - g.interp' * (g.W .* g.ny .* gprS);
  rhspr(vmapDG) = rhspr(vmapDG) + g.interp' * (g.W .* (gvnM - gvnS));

  %%%%%%%%%%%%%%%%%
  % Volume Update %
  %%%%%%%%%%%%%%%%%

  % c.V            interpolates to the cubature nodes
  % c.Dr and c.D.s interpolate  to the cubature nodes and take derivates
  % c.W            is Gauss weights and the Jacobian
  c = OP.dg.cubature;

  cpr = c.V * pr(vmapDG);
  SxTpr = c.Dr' * (c.W .* c.rx .* cpr) + ...
          c.Ds' * (c.W .* c.sx .* cpr);
  SyTpr = c.Dr' * (c.W .* c.ry .* cpr) + ...
          c.Ds' * (c.W .* c.sy .* cpr);

  Sxv1 =  c.V' * ((c.W .* c.rx .* (c.Dr * v1(vmapDG))) + ...
                  (c.W .* c.sx .* (c.Ds * v1(vmapDG))));
  Syv2 =  c.V' * ((c.W .* c.ry .* (c.Dr * v2(vmapDG))) + ...
                  (c.W .* c.sy .* (c.Ds * v2(vmapDG))));

  rhsv1(vmapDG) = rhsv1(vmapDG) + SxTpr;
  rhsv2(vmapDG) = rhsv2(vmapDG) + SyTpr;
  rhspr(vmapDG) = rhspr(vmapDG) - (Sxv1 + Syv2);

  %%%%%%%%%%%%%%%%%%%
  % WADG inverse MJ %
  %%%%%%%%%%%%%%%%%%%
  MI = V * V';
  rhsv1(vmapDG) = MI * (c.V' * (c.WJI .* (c.V * (MI * rhsv1(vmapDG)))));
  rhsv2(vmapDG) = MI * (c.V' * (c.WJI .* (c.V * (MI * rhsv2(vmapDG)))));
  rhspr(vmapDG) = MI * (c.V' * (c.WJI .* (c.V * (MI * rhspr(vmapDG)))));
end
