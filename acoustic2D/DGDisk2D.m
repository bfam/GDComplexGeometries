% NOTE: p should be specified before executing script
r = [2, 2, 2, 5, 5];
sval = [3,3,3,3,3];

format long e
sol_err = [];
for s = 0:sval(p)
  [T, err] = DGDisk2D_func(p, s, 'disk', r(p));
  sol_err(s+1) = err(end);
  if(s > 0)
    disp('DGDisk2D')
    disp('p')
    disp(p)
    disp('error')
    disp(sol_err)
    rate = (log(sol_err(1:end-1)) - log(sol_err(2:end))) / log(2);
    disp('rate')
    disp(rate)
  end
end
disp('FINAL: DGDisk2D')
disp('p')
disp(p)
disp('error')
disp(sol_err)
rate = (log(sol_err(1:end-1)) - log(sol_err(2:end))) / log(2);
disp('rate')
disp(rate)
