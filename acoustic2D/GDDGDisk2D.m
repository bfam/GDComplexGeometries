% NOTE: p should be specified before executing script
r = [2, 2, 2, 5, 5];
N = [9, 14, 19, 24, 28];
sval = [3,3,3,3,3];

format long e
sol_err = [];
for s = 0:sval(p)
  [T, err] = GDDGDisk2D_func(p, N(p), s, 'diskhole', r(p));
  sol_err(s+1) = err(end);
  if(s > 0)
    disp('GDDGDisk2D')
    disp('p')
    disp(p)
    disp('error')
    disp(sol_err)
    rate = (log(sol_err(1:end-1)) - log(sol_err(2:end))) / log(2);
    disp('rate')
    disp(rate)
  end
end
disp('FINAL: GDDGDisk2D')
disp('p')
disp(p)
disp('error')
disp(sol_err)
rate = (log(sol_err(1:end-1)) - log(sol_err(2:end))) / log(2);
disp('rate')
disp(rate)
