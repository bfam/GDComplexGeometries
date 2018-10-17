if exist('p') == 0
  p = 1;
end
if exist('do_plot') == 0
  do_plot = false;
end
if exist('sr') == 0
  sr = 0;
end

Globals2D_gddg
addpath([NODAL_DG_ROOT, '/nodal-dg/Codes1.1/Codes1D'])
addpath([NODAL_DG_ROOT, '/nodal-dg/Codes1.1/Codes2D'])
addpath([NODAL_DG_ROOT, '/nodal-dg/Codes1.1/ServiceRoutines'])
addpath ../src
addpath ../geo/inclusion

Globals2D;


% order of gd approximation

% order of the DG
N = 2*p+1;

% order of the mortar
Nm = 2*p+1;

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

intC = 2*(N+1);
OP.dg.cubature = DGMetricStorageUpdate(CubatureVolumeMesh2D(intC));
%{
OP.dg.cubature = CubatureVolumeMesh2D(intC);
OP.dg.cubature.WJI = OP.dg.cubature.W ./ (OP.dg.cubature.J).^2;
%}

intG = 2*(N+1);
OP.dg.gauss = GaussFaceMesh2D(intG);
[OP.dg.gauss.mapB, OP.dg.gauss.mapMor] = ...
  BuildBCMaps2D_cubature(intC, 5*BCType, OP.dg.gauss.mapM);

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
end
OP.x1 = x1; OP.x2 = x2;

% set up the coupling
% zeros means boundary conditions
% negative face means flipped face
D = 9;
OP.B{1}.toB = [0 2 0 4];
OP.B{1}.toF = [2 1 4 3];

OP.B{2}.toB = [1 3 0 D];
OP.B{2}.toF = [2 1 4 3];

OP.B{3}.toB = [2 0 0 5];
OP.B{3}.toF = [2 1 4 3];

OP.B{4}.toB = [0 D 1 6];
OP.B{4}.toF = [2 1 4 3];

OP.B{5}.toB = [D 0 3 8];
OP.B{5}.toF = [2 1 4 3];

OP.B{6}.toB = [0 7 4 0];
OP.B{6}.toF = [2 1 4 3];

OP.B{7}.toB = [6 8 D 0];
OP.B{7}.toF = [2 1 4 3];

OP.B{8}.toB = [7 0 5 0];
OP.B{8}.toF = [2 1 4 3];

OP.B{9}.toB = [4 5 2 7];
OP.B{9}.toF = [2 1 4 3];
OP.B{9}.isdg = true;
OP.B{9}.quad_order = Nfp-1;
OP.B{9}.vmap = vmapDG;

% Create mortar grids for the DG block

% Here we assume that the "mortar" is aligned with one of the coordinate
% directions with the mapping that the DG faces match the GD faces.
for m = 1:length(vmapMor)
  clear f

  % First we sort the maps so that they match the expected GD ordering
  for k = 1:size(vmapMor{m},2)
    if m == 1 || m == 2
      [~, I] = sort(y(vmapMor{m}(:,k)));
    else
      [~, I] = sort(x(vmapMor{m}(:,k)));
    end
    vmapMor{m}(:,k) = vmapMor{m}(I,k);
    mapMor{m}(:,k) = mapMor{m}(I,k);
  end
  if m == 1 || m == 2
     [~, I] = sort(y(vmapMor{m}(1,:)));
  else
     [~, I] = sort(x(vmapMor{m}(1,:)));
  end
  vmapMor{m} = vmapMor{m}(:,I);
  mapMor{m} = mapMor{m}(:,I);

  % Now we create the data needed to create the mortar
  % (see gd_setup_curved)
  f.x1 = x(vmapMor{m}(:));
  f.x2 = y(vmapMor{m}(:));
  f.corners = [f.x1(  1,   1), f.x2(  1,   1);
               f.x1(end, end), f.x2(end, end)];
  if m == 1 || m == 2
     f.g  = [y(vmapMor{m}(1,:)),y(vmapMor{m}(end,end))]';
     f.rq = f.x2(:);
  else
     f.g = [x(vmapMor{m}(1,:)),x(vmapMor{m}(end,end))]';
     f.rq = f.x1(:);
  end
  f.rq = 2 * (f.rq - f.g(1)) / (f.g(end) - f.g(1)) - 1;
  f.g  = 2 * (f.g  - f.g(1)) / (f.g(end) - f.g(1)) - 1;
  f.P = 1;
  f.vmap = vmapMor{m}(:);

  OP.B{9}.f{m} = f;
end

% OP.dg.gauss.mapB = sort(unique([OP.dg.gauss.mapB;OP.dg.gauss.mapMor]));
% OP.dg.gauss.mapMor = [];

OP.B = mortar_create(OP.B);
if do_plot
  plot_mesh(OP.B)
  hold on
  PlotMesh2D
  hold off
end

cx = cos(linspace(0, 2*pi, 101));
cy = sin(linspace(0, 2*pi, 101));
c2_x = cyl.r(1) * cx + cyl.x(1);
c2_y = cyl.r(1) * cy + cyl.y(1);
c3_x = cyl.r(2) * cx + cyl.x(2);
c3_y = cyl.r(2) * cy + cyl.y(2);
c4_x = cyl.r(3) * cx + cyl.x(3);
c4_y = cyl.r(3) * cy + cyl.y(3);
c5_x = cyl.r(4) * cx + cyl.x(4);
c5_y = cyl.r(4) * cy + cyl.y(4);

Exact2D = @(t,x,y) deal(zeros(size(x)), zeros(size(y)),  ...
                        exp(-2*(((hypot(x,y)-10).^2))));
[v1, v2, pr] = L2ProjExact(OP, 0);
if do_plot
  subplot(2,3,1)
  gd_contourf(OP.B, pr, false, 100);
  cmap
  caxis([-1 1])
  hold on
  PlotField2D(N, x, y, pr(vmapDG));
  plot_mesh(OP.B)
  hold on
  PlotMesh2D
  plot(c2_x, c2_y, 'k')
  plot(c3_x, c3_y, 'k')
  plot(c4_x, c4_y, 'k')
  plot(c5_x, c5_y, 'k')
  hold off
  axis image
  drawnow
end

clear Exact2D
dt = compute_dt(OP)/4;
StepTime = 1; Ntsteps = ceil(StepTime/dt); dt = StepTime/Ntsteps;
for k = 1:30
  [v1,v2,pr,T] = Acoustic2D(OP, v1, v2, pr, dt, Ntsteps, 2*p+2);
  if do_plot
    subplot(2,3,mod(k, 6)+1)
    gd_contourf(OP.B, pr, false, 100);
    cmap
    caxis([-1 1])
    hold on
    PlotField2D(N, x, y, pr(vmapDG));
    % PlotMesh2D()
    plot(c2_x, c2_y, 'k')
    plot(c3_x, c3_y, 'k')
    plot(c4_x, c4_y, 'k')
    plot(c5_x, c5_y, 'k')
    hold off
    axis image
    colorbar
    title( sprintf('n%d :: r%d :: t%d', 2*p+1, sr, k))
    drawnow
  end
  eval(sprintf('save data_n%d_r%d_t%d v1 v2 pr x1 x2', 2*p+1, sr, k))
end

