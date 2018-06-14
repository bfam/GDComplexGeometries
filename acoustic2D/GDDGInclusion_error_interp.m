clear

Globals2D_gddg
addpath([NODAL_DG_ROOT, '/nodal-dg/Codes1.1/Codes1D'])
addpath([NODAL_DG_ROOT, '/nodal-dg/Codes1.1/Codes2D'])
addpath([NODAL_DG_ROOT, '/nodal-dg/Codes1.1/ServiceRoutines'])
addpath ../src
addpath ../geo/inclusion

Globals2D;

% order of gd approximation
p = 2;

% order of the mortar
Nm = 2*p+1;

% order of the DG
N = 2*p+1;

[r_q, ~] = gauss_cofs(Nm+2);
r_q2 = [r_q-1;r_q+1]/2;
r_q4 = [r_q2-1;r_q2+1]/2;
Iq_q2 = lagrange_interpolation_matrix(r_q, r_q2);
Iq_q4 = lagrange_interpolation_matrix(r_q, r_q4);

for sr = 2:-1:0
  [Nv, VX, VY, K, EToV, BCType] = MeshReaderGmshBC2D('inclusion_r0.msh');
  MakeInclusionsVerts2D()
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

  R = 5;
  L = 15;

  h0 = 0.75;
  N2R = 8*N;
  h0 = (2*R / N2R);
  NLR = ceil((L - R) / h0);

  N2R = (2^sr) * N2R;
  NLR = (2^sr) * NLR;

  grid{1}.x1 = [-L -R];
  grid{1}.x2 = [-L -R];
  grid{1}.N1 = NLR; grid{1}.N2 = NLR;

  grid{2}.x1 = [-R R];
  grid{2}.x2 = [-L -R];
  grid{2}.N1 = N2R; grid{2}.N2 = NLR;

  grid{3}.x1 = [R L];
  grid{3}.x2 = [-L -R];
  grid{3}.N1 = NLR; grid{3}.N2 = NLR;

  grid{4}.x1 = [-L -R];
  grid{4}.x2 = [-R R];
  grid{4}.N1 = NLR; grid{4}.N2 = N2R;

  grid{5}.x1 = [R L];
  grid{5}.x2 = [-R R];
  grid{5}.N1 = NLR; grid{5}.N2 = N2R;

  grid{6}.x1 = [-L -R];
  grid{6}.x2 = [R L];
  grid{6}.N1 = NLR; grid{6}.N2 = NLR;

  grid{7}.x1 = [-R R];
  grid{7}.x2 = [R L];
  grid{7}.N1 = N2R; grid{7}.N2 = NLR;

  grid{8}.x1 = [R L];
  grid{8}.x2 = [R L];
  grid{8}.N1 = NLR; grid{8}.N2 = NLR;

  %
  % Initialize the blocks
  %

  eval(sprintf('load data/data_n%d_r%d_t15', 2*p+1, sr))
  % eval(sprintf('load data_n%d_r%d_t0', 2*p+1, sr))
  xc1 = x1;
  xc2 = x2;

  x1 = x(:); x2 = y(:);

  for k = 1:length(grid)
    % OP.B{k} = gd_setup_curved(grid{k}.N1, grid{k}.N2, p, Nm+1, grid{k});
    OP.B{k} = gd_setup_affine(grid{k}.N1, grid{k}.N2, p, Nm+1, ...
                              grid{k}.x1, grid{k}.x2, 0);

    % Set up the global vmap
    OP.B{k}.vmap = length(x1) + (1:length(OP.B{k}.x1))';
    OP.B{k}.vmap = reshape(OP.B{k}.vmap, OP.B{k}.Np2, OP.B{k}.Np1);

    % Set up the global grid
    x1 = [x1;OP.B{k}.x1]; x2 = [x2;OP.B{k}.x2];

    x1_f{sr+1,k} = reshape(OP.B{k}.x1, OP.B{k}.N2+1, OP.B{k}.N1+1);
    x1_f{sr+1,k} = OP.B{k}.P2 * (OP.B{k}.E2 * x1_f{sr+1,k});
    x1_f{sr+1,k} = (x1_f{sr+1,k} * OP.B{k}.E1') * OP.B{k}.P1';

    x2_f{sr+1,k} = reshape(OP.B{k}.x2, OP.B{k}.N2+1, OP.B{k}.N1+1);
    x2_f{sr+1,k} = OP.B{k}.P2 * (OP.B{k}.E2 * x2_f{sr+1,k});
    x2_f{sr+1,k} = (x2_f{sr+1,k} * OP.B{k}.E1') * OP.B{k}.P1';

    pr_f{sr+1,k} = pr(OP.B{k}.vmap);
    pr_f{sr+1,k} = OP.B{k}.P2 * (OP.B{k}.E2 * pr_f{sr+1,k});
    pr_f{sr+1,k} = (pr_f{sr+1,k} * OP.B{k}.E1') * OP.B{k}.P1';

    v1_f{sr+1,k} = v1(OP.B{k}.vmap);
    v1_f{sr+1,k} = OP.B{k}.P2 * (OP.B{k}.E2 * v1_f{sr+1,k});
    v1_f{sr+1,k} = (v1_f{sr+1,k} * OP.B{k}.E1') * OP.B{k}.P1';

    v2_f{sr+1,k} = v2(OP.B{k}.vmap);
    v2_f{sr+1,k} = OP.B{k}.P2 * (OP.B{k}.E2 * v2_f{sr+1,k});
    v2_f{sr+1,k} = (v2_f{sr+1,k} * OP.B{k}.E1') * OP.B{k}.P1';
    if sr == 0
      x1_f{sr+1,k} = kron(speye(OP.B{k}.N2), Iq_q4) * x1_f{sr+1,k};
      x1_f{sr+1,k} = x1_f{sr+1,k} * kron(speye(OP.B{k}.N1), Iq_q4');

      x2_f{sr+1,k} = kron(speye(OP.B{k}.N2), Iq_q4) * x2_f{sr+1,k};
      x2_f{sr+1,k} = x2_f{sr+1,k} * kron(speye(OP.B{k}.N1), Iq_q4');

      pr_f{sr+1,k} = kron(speye(OP.B{k}.N2), Iq_q4) * pr_f{sr+1,k};
      pr_f{sr+1,k} = pr_f{sr+1,k} * kron(speye(OP.B{k}.N1), Iq_q4');

      v1_f{sr+1,k} = kron(speye(OP.B{k}.N2), Iq_q4) * v1_f{sr+1,k};
      v1_f{sr+1,k} = v1_f{sr+1,k} * kron(speye(OP.B{k}.N1), Iq_q4');

      v2_f{sr+1,k} = kron(speye(OP.B{k}.N2), Iq_q4) * v2_f{sr+1,k};
      v2_f{sr+1,k} = v2_f{sr+1,k} * kron(speye(OP.B{k}.N1), Iq_q4');
    elseif sr == 1
      x1_f{sr+1,k} = kron(speye(OP.B{k}.N2), Iq_q2) * x1_f{sr+1,k};
      x1_f{sr+1,k} = x1_f{sr+1,k} * kron(speye(OP.B{k}.N1), Iq_q2');

      x2_f{sr+1,k} = kron(speye(OP.B{k}.N2), Iq_q2) * x2_f{sr+1,k};
      x2_f{sr+1,k} = x2_f{sr+1,k} * kron(speye(OP.B{k}.N1), Iq_q2');

      pr_f{sr+1,k} = kron(speye(OP.B{k}.N2), Iq_q2) * pr_f{sr+1,k};
      pr_f{sr+1,k} = pr_f{sr+1,k} * kron(speye(OP.B{k}.N1), Iq_q2');

      v1_f{sr+1,k} = kron(speye(OP.B{k}.N2), Iq_q2) * v1_f{sr+1,k};
      v1_f{sr+1,k} = v1_f{sr+1,k} * kron(speye(OP.B{k}.N1), Iq_q2');

      v2_f{sr+1,k} = kron(speye(OP.B{k}.N2), Iq_q2) * v2_f{sr+1,k};
      v2_f{sr+1,k} = v2_f{sr+1,k} * kron(speye(OP.B{k}.N1), Iq_q2');
    end
  end
  if norm(x1-xc1,'inf') || norm(x2-xc2,'inf')
    warning(sprintf(['possible problem with grid: '...
                     '||x1-xc1||_{inf} = %e, ||x2-xc2||_{inf} = %e'], ...
                     norm(x1-xc1,'inf'), norm(x2-xc2,'inf')))
  end
  pr_dg{sr+1} = RefineUniformField2D(2-sr, pr(vmapDG));
  v1_dg{sr+1} = RefineUniformField2D(2-sr, v1(vmapDG));
  v2_dg{sr+1} = RefineUniformField2D(2-sr, v2(vmapDG));
  x1_dg{sr+1} = RefineUniformField2D(2-sr, x1(vmapDG));
  x2_dg{sr+1} = RefineUniformField2D(2-sr, x2(vmapDG));

  % if not the finest mesh, need to correctly interpolate the curved elements
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
          end
        end
        I = Vandermonde2D(N, rf, sf) * invV;

        pr_dg{sr+1}(:,(k-1)*ne + (1:ne)) = reshape(I * pr(vmapDG(:,k)), Np, ne);
        v1_dg{sr+1}(:,(k-1)*ne + (1:ne)) = reshape(I * v1(vmapDG(:,k)), Np, ne);
        v2_dg{sr+1}(:,(k-1)*ne + (1:ne)) = reshape(I * v2(vmapDG(:,k)), Np, ne);
        x1_dg{sr+1}(:,(k-1)*ne + (1:ne)) = reshape(I * x(:,k), Np, ne);
        x2_dg{sr+1}(:,(k-1)*ne + (1:ne)) = reshape(I * y(:,k), Np, ne);
      end
    end
  else
    intC = 2*(N+1);
    cubature = CubatureVolumeMesh2D(intC);
    cubature.WJI = cubature.W ./ (cubature.J).^2;
    for k = 1:length(OP.B)
      WJ{k} = OP.B{k}.w .* OP.B{k}.J;
    end
  end
end

err_gd = [0,0];
err_dg = [0,0];
for sr = 0:1
  for k = 1:length(WJ)
    dv1 = v1_f{sr+1,k}(:) - v1_f{sr+2,k}(:);
    dv2 = v2_f{sr+1,k}(:) - v2_f{sr+2,k}(:);
    dpr = pr_f{sr+1,k}(:) - pr_f{sr+2,k}(:);
    err_gd(sr+1) = err_gd(sr+1) + sum(WJ{k} .* dv1.^2) / 2;
    err_gd(sr+1) = err_gd(sr+1) + sum(WJ{k} .* dv2.^2) / 2;
    err_gd(sr+1) = err_gd(sr+1) + sum(WJ{k} .* dpr.^2) / 2;
  end

  % Add in DG portion
  c = cubature;
  dv1 = c.V*v1_dg{sr+1} - c.V*v1_dg{sr+2};
  dv2 = c.V*v2_dg{sr+1} - c.V*v2_dg{sr+2};
  dpr = c.V*pr_dg{sr+1} - c.V*pr_dg{sr+2};
  err_dg(sr+1) = err_dg(sr+1) + sum(sum(diag(c.w) * (c.J .* dv1.^2))) / 2;
  err_dg(sr+1) = err_dg(sr+1) + sum(sum(diag(c.w) * (c.J .* dv2.^2))) / 2;
  err_dg(sr+1) = err_dg(sr+1) + sum(sum(diag(c.w) * (c.J .* dpr.^2))) / 2;

  dx1 = c.V*x1_dg{sr+1} - c.V*x1_dg{sr+2};
  err_dg_x1(sr+1) = sum(sum(diag(c.w) * (c.J .* dx1.^2)));
  dx2 = c.V*x2_dg{sr+1} - c.V*x2_dg{sr+2};
  err_dg_x2(sr+1) = sum(sum(diag(c.w) * (c.J .* dx2.^2)));
end

err = sqrt(err_gd + err_dg);
err_gd = sqrt(err_gd);
err_dg = sqrt(err_dg);

disp(log2(err(1)    / err(2)   ))
disp(log2(err_gd(1) / err_gd(2)))
disp(log2(err_dg(1) / err_dg(2)))

err_dg_x1 = sqrt(err_dg_x1);
err_dg_x2 = sqrt(err_dg_x2);
disp(log2(err_dg_x1(1) / err_dg_x1(2)))
disp(log2(err_dg_x2(1) / err_dg_x2(2)))
