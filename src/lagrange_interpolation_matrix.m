% Compute the Lagrange interpolation matrix from r_src to r_dst (uses
% barycentric interpolation formula)
function I = lagrange_interpolation_matrix(r_src, r_dst)
  w_src = lagrange_weights(r_src);

  I = zeros(length(r_dst), length(r_src));

  for k = 1:length(r_dst)
    for j = 1:length(r_src)
      I(k, j) = w_src(j) / (r_dst(k) - r_src(j));
      if(~isfinite(I(k,j)))
        I(k, :) = 0;
        I(k, j) = 1;
        break
      end
    end
    d = sum(I(k, :));
    I(k, :) = I(k, :) / d;
  end

end
