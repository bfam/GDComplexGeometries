% Create the mortar for GD operator stored as cell array B
% Modifies the face struct of connected blocks
function B = mortar_create(B)
  % Check and fix connectivity:
  for km = 1:length(B)
    for fm = 1:4
      kp = B{km}.toB(fm);
      fp = B{km}.toF(fm);
      if(kp > 0)
        % Make sure the sides match
        if km ~= B{kp}.toB(fp) || sign(fp)*fm ~= B{kp}.toF(fp)
          error(sprintf('problem with block %d and face %d', km, fm))
        end

        % First check the conforming faces to make sure they are OK
        % TODO: Fix this sort that we check the points too (not just length)
        keep_working = true;
        if length(B{km}.f{fm}.x1) == length(B{kp}.f{fp}.x1)
          keep_working = false;
          mx = 0;
          mx = max(mx,norm(B{km}.f{fm}.x1 - B{kp}.f{fp}.x1));
          mx = max(mx,norm(B{km}.f{fm}.x2 - B{kp}.f{fp}.x2));
          if mx > sqrt(eps)
            % Try to flip
            mx = 0;
            mx = max(mx,norm(B{km}.f{fm}.x1 - B{kp}.f{fp}.x1(end:-1:1)));
            mx = max(mx,norm(B{km}.f{fm}.x2 - B{kp}.f{fp}.x2(end:-1:1)));
            if mx < sqrt(eps)
              flds = fieldnames(B{kp}.f{fp});
              for k = 1:length(flds)
                if ~strcmp('vmap', flds{k})
                  eval(sprintf('B{kp}.f{fp}.%s = flipud(B{kp}.f{fp}.%s);',...
                    flds{k}, flds{k}));
                end
              end
            else
              keep_working = true;
            end
          end
        end
        if keep_working
          % Now handle the non-conforming faces

          % First check the orientation and flip if needed
          shift = [0,0];
          if isfield(B{km}, 'toShftX')
            shift = [B{km}.toShftX(fm), B{km}.toShftY(fm)];
          end
          mx = norm(B{km}.f{fm}.corners - B{kp}.f{fp}.corners + [shift;shift]);
          if mx > sqrt(eps)
            flds = fieldnames(B{kp}.f{fp});

            % We don't flip x1 and x2 because we project them below
            for k = 1:length(flds)
              if ~strcmp('vmap', flds{k}) && ...
                 ~strcmp('x1', flds{k}) && ~strcmp('x2', flds{k})
                eval(sprintf('B{kp}.f{fp}.%s = flipud(B{kp}.f{fp}.%s);',...
                  flds{k}, flds{k}));
              end
            end

            mx = norm(B{km}.f{fm}.corners - B{kp}.f{fp}.corners + [shift;shift]);
            if mx > sqrt(eps)
              error(sprintf('problem with block %d and face %d', km, fm))
            end

            g = B{kp}.f{fp}.g;
            B{kp}.f{fp}.g = g(1) - g + g(end);
            rq = B{kp}.f{fp}.rq;
            B{kp}.f{fp}.rq = rq(1) - rq + rq(end);
          end
          m.g  = B{km}.f{fm}.g;
          m.r  = B{km}.f{fm}.rq;
          m.Pg = B{km}.f{fm}.P;

          p.g  = B{kp}.f{fp}.g;
          p.r  = B{kp}.f{fp}.rq;
          p.Pg = B{kp}.f{fp}.P;
          g = mortar_setup(p, m, max(B{km}.quad_order, B{kp}.quad_order)+1);

          B{km}.f{fm}.rq = g.r(:);
          B{km}.f{fm}.w = g.w(:);
          B{km}.f{fm}.P = g.Im;

          B{kp}.f{fp}.rq = g.r(:);
          B{kp}.f{fp}.w  = g.w(:);
          B{kp}.f{fp}.P  = g.Ip;

          % TODO: Should we do something besides average?
          x1 = shift(1) + (B{km}.f{fm}.P*B{km}.f{fm}.x1 + B{kp}.f{fp}.P*B{kp}.f{fp}.x1)/2;
          x2 = shift(2) + (B{km}.f{fm}.P*B{km}.f{fm}.x2 + B{kp}.f{fp}.P*B{kp}.f{fp}.x2)/2;
          %{
          clf;
          subplot(2,1,1)
          plot(B{km}.f{fm}.P*B{km}.f{fm}.x1 - B{kp}.f{fp}.P*B{kp}.f{fp}.x1, '*')
          hold on
          plot(B{km}.f{fm}.P*B{km}.f{fm}.x2 - B{kp}.f{fp}.P*B{kp}.f{fp}.x2, '*')
          hold off
          subplot(2,1,2)
          plot(x1)
          hold on
          plot(x2)
          hold off
          title([km, fm])
          pause
          %}
          Dx1 = g.D * x1; Dx2 = g.D * x2;
          sJ = sqrt(Dx1.^2 + Dx2.^2);

          B{km}.f{fm}.x1 = x1;
          B{km}.f{fm}.x2 = x2;
          B{km}.f{fm}.sJ = sJ;
          if fm == 1 || fm == 4
            B{km}.f{fm}.n1 = -Dx2 ./ sJ;
            B{km}.f{fm}.n2 =  Dx1 ./ sJ;
          elseif fm == 2 || fm == 3
            B{km}.f{fm}.n1 =  Dx2 ./ sJ;
            B{km}.f{fm}.n2 = -Dx1 ./ sJ;
          else
            error(sprintf('problem with block %d and face %d', km, fm))
          end
          B{kp}.f{fp}.sJ =  B{km}.f{fm}.sJ;
          B{kp}.f{fp}.n1 = -B{km}.f{fm}.n1;
          B{kp}.f{fp}.n2 = -B{km}.f{fm}.n2;
          B{kp}.f{fp}.x1 =  B{km}.f{fm}.x1;
          B{kp}.f{fp}.x2 =  B{km}.f{fm}.x2;
        end
      end
    end
  end
end

function g = mortar_setup(p, m, Npm)
  if(~isfield(p, 'Pg'))
    p.Pg = 1;
  end
  if(~isfield(m, 'Pg'))
    m.Pg = 1;
  end

  if( ~issorted(p.g) || ~issorted(p.r) || ...
      ~issorted(m.g) || ~issorted(m.r) )
     error('something on the mortar not properly sorted')
  end

  [g.g, g.r, g.w, g.D] = mortar_combine_grid(p.g, m.g, Npm);
  g.Im = mortar_interpolation(m.r, m.g, g.r, g.g, m.Pg);
  g.Ip = mortar_interpolation(p.r, p.g, g.r, g.g, p.Pg);
end

function [g, r, w, D] = mortar_combine_grid(gp, gm, Nm)
  % g  end points of mortar cells
  % r  Gauss nodes   for mortar grid (order Nm - 1)
  % w  Gauss weights for mortar grid (order Nm - 1)

  r = 0;
  w = 0;

  g = sort([gp;gm]);
  tol = sqrt(eps);
  if(abs(gp(1)-gm(1)) > tol || abs(gp(end)-gm(end)) > tol)
    error('grids must match at end points')
  end
  g = g([diff(g)>tol; true]); % true at end needed so we don't miss last point

  [rgl, wgl] = gauss_cofs(Nm);
  Dgl = spectral_derivative(rgl);
  D = kron(2 * diag(sparse(1./diff(g'))), Dgl);

  r = 0.5*(rgl+1)*diff(g') + ones(size(rgl))*g(1:end-1)';
  w = 0.5*wgl*diff(g');

  %{
  plot(r,ones(size(r,1),1)*max(r),'.')
  hold on
  plot(g(1:end-1),g(2:end),'k*')
  plot(g(2:end  ),g(2:end),'k*')
  hold off
  %}
end

function I = mortar_interpolation(r_src, g_src, r_dst, g_dst, Pg)

  if(size(r_src,1) == 1 || size(r_src,2)==1)
    Np_src = length(r_src) / (length(g_src)-1);
    r_src  = reshape(r_src, Np_src, length(g_src)-1);
  end

  if(size(r_dst,1) == 1 || size(r_dst,2)==1)
    Np_dst = length(r_dst) / (length(g_dst)-1);
    r_dst  = reshape(r_dst, Np_dst, length(g_dst)-1);
  end

  Np_src = size(r_src,1);
  Np_dst = size(r_dst,1);

  if(Np_src > Np_dst)
    warning('destination space seems to be lower order!')
  end

  I_I = zeros(Np_dst, Np_src, length(g_dst)-1);
  I_J = zeros(Np_dst, Np_src, length(g_dst)-1);
  Ikm = zeros(Np_dst, Np_src, length(g_dst)-1);

  I_I(:,:,1) = (1:Np_dst)'*ones(1,Np_src);
  I_J(:,:,1) = ones(Np_dst,1)*(1:Np_src);

  k = 1;
  tol = sqrt(eps);
  for m = 1:length(g_dst)-1
    % Find the overlap
    while(g_src(k+1) <= g_dst(m))
       k = k+1;
    end
    if((g_src(k) - g_dst(m)) > tol || (g_dst(m+1) - g_src(k+1)) > tol)
      error('destination not nested in source')
    end

    % store the indices and coefficients for the interpolation
    I_I(:,:,m) = Np_dst*(m-1) + I_I(:,:,1);
    I_J(:,:,m) = Np_src*(k-1) + I_J(:,:,1);
    Ikm(:,:,m) = lagrange_interpolation_matrix(r_src(:,k), r_dst(:,m));
  end

  % build the interpolation matrix
  I = sparse(I_I(:),I_J(:),Ikm(:));

  % norm(I * r_src(:) - r_dst(:))

  I = I * Pg;

end
