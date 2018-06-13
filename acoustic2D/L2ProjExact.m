% Approximate L2 projection of exact solution at time time
function [v1, v2, pr] = L2ProjExact(OP, time)

Globals2D_gddg;

v1 = zeros(size(OP.x1));
v2 = zeros(size(OP.x1));
pr = zeros(size(OP.x1));


%%%%%%%%%%%%%%%%%%%%%%%
% Galerkin Difference %
%%%%%%%%%%%%%%%%%%%%%%%

err_gd = 0;
eng_gd = 0;

for k = 1:length(OP.B)

  if isfield(OP.B{k}, 'isdg')
    continue
  end

  x1g = OP.B{k}.P1F*(OP.B{k}.P2F*(OP.B{k}.x1));
  x2g = OP.B{k}.P1F*(OP.B{k}.P2F*(OP.B{k}.x2));

  [ev1, ev2, epr] = Exact2D(time, x1g, x2g);

  MJev1 = reshape((OP.B{k}.P2F' * (OP.B{k}.P1F' * (OP.B{k}.w .* OP.B{k}.J .* ev1))), size(OP.B{k}.vmap));
  MJev2 = reshape((OP.B{k}.P2F' * (OP.B{k}.P1F' * (OP.B{k}.w .* OP.B{k}.J .* ev2))), size(OP.B{k}.vmap));
  MJepr = reshape((OP.B{k}.P2F' * (OP.B{k}.P1F' * (OP.B{k}.w .* OP.B{k}.J .* epr))), size(OP.B{k}.vmap));

  v1(OP.B{k}.vmap) = OP.B{k}.massInv(MJev1);
  v2(OP.B{k}.vmap) = OP.B{k}.massInv(MJev2);
  pr(OP.B{k}.vmap) = OP.B{k}.massInv(MJepr);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discontinuous Galerkin %
%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(OP, 'dg')
  Globals2D;

  c = OP.dg.cubature;

  x1g = c.V*OP.x1(vmapDG);
  x2g = c.V*OP.x2(vmapDG);

  [ev1, ev2, epr] = Exact2D(time, x1g, x2g);

  MJev1 = c.V' * (c.W .* ev1);
  MJev2 = c.V' * (c.W .* ev2);
  MJepr = c.V' * (c.W .* epr);



  MI = V * V';
  v1(vmapDG) = MI * (c.V' * (c.WJI .* (c.V * (MI * MJev1))));
  v2(vmapDG) = MI * (c.V' * (c.WJI .* (c.V * (MI * MJev2))));
  pr(vmapDG) = MI * (c.V' * (c.WJI .* (c.V * (MI * MJepr))));

end
