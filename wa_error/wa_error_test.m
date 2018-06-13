clear
addpath ../src
format short e

% Starting figure number
fig = 0;

% set to true to do plotting
do_plotting = false;

% set to true dump tikz data
do_tikz = false;

% beta value for the underlying coordinate transform
for beta = [0.22, 1/8]
  % Ng = 0 is extrapolation method and Ng = 1 is ghost basis
  for Ng = [0, 1]
    % mode of the function being approximated
    for k = [1, 10] * pi /2
      % GD half order (actual GD polynomial order is 2*p+1
      for p = 2

        fprintf(['[beta, Ng, k, p] = [', num2str([beta, Ng, k, p]), ']\n'])
        disp('        level     error L2     error WA    error ~WA')

        % Number of quadrature points to use for "exact" integration
        Np_q = 2 * p + 2;

        % N+1 is number of interior GD points
        N = (2*p+1)*2.^(0:3);
        for j = 1:length(N)
          [errL2(j), errWA(j), errJnWA(j)] = wa_error_fun(beta, N(j), N(j), p, Ng*p, Np_q, k);
          disp([j errL2(j), errWA(j), errJnWA(j)])
        end

        if do_plotting
          fig = fig + 1;
          figure(fig)
          h = 1:length(N);
          semilogy(h, errL2, '*-', h, errWA, '*-', h, errJnWA, '*-', ...
          h, errWA-errL2, '*--', h, errJnWA-errL2, '*--')
          legend('L2', 'WA', 'JnWA', 'WA-L2', 'JnWA-L2')
          axis tight
          title(['[beta, Ng, k, p] = [', num2str([beta, Ng, k, p]), ']'])
          drawnow
        end

        if do_tikz
          fprintf('        %% beta = %e\n', beta)
          fprintf('        %% mode k = %d\n', k);
          fprintf('        %% poly order = %d\n', 2*p+1)
          fprintf('        %% num ghost = %d\n', Ng*p)
          fprintf('        %% num gauss quad pnts = %d\n', Np_q)
          fprintf('        \\addplot[color=grn,mark=triangle*,line width=1] plot coordinates {%%\n')
          for j = 1:length(N)
            fprintf('          (%e, %e)\n', h(j), errL2(j))
          end
          fprintf('        };\n')
          fprintf('        \\addplot[color=red,mark=diamond*,line width=1] plot coordinates {%%\n')
          for j = 1:length(N)
            fprintf('          (%e, %e)\n', h(j), errWA(j))
          end
          fprintf('        };\n')
          fprintf('        \\addplot[color=grn,mark=square*,line width=1] plot coordinates {%%\n')
          for j = 1:length(N)
            fprintf('          (%e, %e)\n', h(j), errJnWA(j))
          end
          fprintf('        };\n')
            fprintf('        \\addplot[color=red,mark=pentagon*,line width=1,dashed] plot coordinates {%%\n')
          for j = 1:length(N)
            fprintf('          (%e, %e)\n', h(j), errWA(j)-errL2(j))
          end
          fprintf('        };\n')
          fprintf('        \\addplot[color=blu,mark=*,line width=1,dashed] plot coordinates {%%\n')
          for j = 1:length(N)
            fprintf('          (%e, %e)\n', h(j), errJnWA(j)-errL2(j))
          end
          fprintf('        };\n\n\n')
        end
      end
    end
  end
end

if do_tikz
  fprintf('        \\legend{L2-projection\\\\WA-projection\\\\%%\n')
  fprintf('                approximate WA-projection\\\\WA error\\\\approximate WA error\\\\}\n')
end
