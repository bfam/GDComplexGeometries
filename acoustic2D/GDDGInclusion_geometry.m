clear

addpath ../nodal-dg/Codes1.1/Codes1D
addpath ../nodal-dg/Codes1.1/Codes2D
addpath ../nodal-dg/Codes1.1/ServiceRoutines
addpath ../src
addpath ../geo/inclusion

Globals2D_gddg
Globals2D;

% order of gd approximation
p = 3;

% order of the mortar
Nm = 2*p+1;

% order of the DG
N = 2*p+1;

for sr = 2:-1:0
  [Nv, VX, VY, K, EToV, BCType] = MeshReaderGmshBC2D('inclusion_r0.msh');
  RefineUniform2D(sr, @MakeInclusionsVerts2D);
  StartUp2D;
  % BuildBCMaps2D

  % setup the curved grid
  cyl.bc = [5, 6, 7, 8];
  cyl.r  = [2, 2.5, 1.5, 1];
  cyl.x  = [2, -2, -3, 2.5];
  cyl.y  = [2, -2, 2.5, -3];
  for c = 1:length(cyl.bc)
    [k, f] = find(BCType == cyl.bc(c));
    if(~isempty(k))
      outfaces = [k,f];
      MakeCylinder2D(outfaces, cyl.r(c), cyl.x(c), cyl.y(c));
    end
    [k, f] = find(BCType == cyl.bc(c));
  end

  BuildBCMaps2D_gddg(cyl.bc);
  vmapDG = reshape([1:Np*K], size(x));

  x1_dg{sr+1} = RefineUniformField2D(2-sr, x);
  x2_dg{sr+1} = RefineUniformField2D(2-sr, y);

  if(sr < 2)

    ne = 4^(2-sr);
    [r0f, s0f] = RefineUniformRS2D(r,s);

    % if we
    if sr == 0
      I = Vandermonde2D(N, r0f, s0f) * invV;
      r0f = I * reshape(r0f, Np, 4); r0f = r0f(:);
      s0f = I * reshape(s0f, Np, 4); s0f = s0f(:);
    end

    for k = 1:size(x,2)
      rf = r0f;
      sf = s0f;
      if ~isempty(intersect(BCType(k,:), cyl.bc))
        for e = 1:ne
          for d = 1:length(r)
            xf = x1_dg{3}(d, (k-1)*ne + e);
            yf = x2_dg{3}(d, (k-1)*ne + e);
            r0 = rf(d + (e-1)*length(r));
            s0 = sf(d + (e-1)*length(r));
            [r1, s1] = PhysicalToReference(xf, yf, r0, s0, x(:, k), y(:, k), ...
                     Dr * x(:, k), Dr * y(:, k), Ds * x(:, k), Ds * y(:, k), ...
                     invV, N);
            rf(d + (e-1)*length(r)) = r1;
            sf(d + (e-1)*length(r)) = s1;
            I = Vandermonde2D(N, r1, s1) * invV;
          end
        end
      end
      I = Vandermonde2D(N, rf, sf) * invV;

      x1_dg{sr+1}(:,(k-1)*ne + (1:ne)) = reshape(I * x(:,k), Np, ne);
      x2_dg{sr+1}(:,(k-1)*ne + (1:ne)) = reshape(I * y(:,k), Np, ne);
      if(sr == 0)
        max(max(abs(x1_dg{sr+1}(:,(k-1)*ne + (1:ne)) - x1_dg{3}(:,(k-1)*ne + (1:ne)))))
        max(max(abs(x2_dg{sr+1}(:,(k-1)*ne + (1:ne)) - x2_dg{3}(:,(k-1)*ne + (1:ne)))))
      end
    end
  else
    intC = 2*(N+1);
    cubature = CubatureVolumeMesh2D(intC);
    cubature.WJI = cubature.W ./ (cubature.J).^2;
  end
end

for sr = 0:1
  % Add in DG portion
  c = cubature;
  dx1 = c.V*x1_dg{sr+1} - c.V*x1_dg{sr+2};
  err_dg_x1(sr+1) = sum(sum(diag(c.w) * (c.J .* dx1.^2)));
  dx2 = c.V*x2_dg{sr+1} - c.V*x2_dg{sr+2};
  err_dg_x2(sr+1) = sum(sum(diag(c.w) * (c.J .* dx2.^2)));
end

err_dg_x1 = sqrt(err_dg_x1);
err_dg_x2 = sqrt(err_dg_x2);
disp(log2(err_dg_x1(1) / err_dg_x1(2)))
disp(log2(err_dg_x2(1) / err_dg_x2(2)))
