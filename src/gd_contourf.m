% Make contour plot
% Inputs:
%   B    :: cell array of GD operators
%   fld  :: data to plot
%   mesh :: (optional) plot the mesh (default: true)
%   N    :: (optional) number of contour lines (default: 20)
function gd_contourf(B,fld,mesh,N)
if nargin < 3
  mesh = true;
end
if nargin < 4
  N = 20;
end

for k = 1:length(B)
  if isfield(B{k},'isdg')
    continue
  end

  x1  = reshape(B{k}.x1, B{k}.Np2, B{k}.Np1);
  x2  = reshape(B{k}.x2, B{k}.Np2, B{k}.Np1);
  fl = reshape(fld(B{k}.vmap), B{k}.Np2, B{k}.Np1);
  p = B{k}.p;

  contourf(B{k}.P2_ge * x1 * B{k}.P1_ge', B{k}.P2_ge * x2 * B{k}.P1_ge', ...
           B{k}.P2_ge * fl * B{k}.P1_ge', N, 'LineStyle', 'none')
  % contourf(x1(p+1:end-p,p+1:end-p), x2(p+1:end-p,p+1:end-p), ...
  %          fl(p+1:end-p,p+1:end-p), N, 'LineStyle', 'none')
  hold on

  if mesh
    if isfield(B{k}, 'E1')
      x1 = B{k}.E2 * x1 * B{k}.E1';
      x2 = B{k}.E2 * x2 * B{k}.E1';
    end
    plot(x1(p+1:end-p,p+1:end-p) , x2(p+1:end-p,p+1:end-p) , 'k');
    plot(x1(p+1:end-p,p+1:end-p)', x2(p+1:end-p,p+1:end-p)', 'k');
  end
end

hold off
axis image
