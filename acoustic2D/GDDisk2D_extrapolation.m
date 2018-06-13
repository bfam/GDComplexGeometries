% NOTE: p should be specified before executing script
r = [2, 2, 2, 5, 5];
N = [4, 7, 10, 14, 16];
sval = [3,3,3,3,3];

format long e
sol_err = [];
for s = 0:sval(p)
  [T, err] = GDDiskDriver2D_func(p, N(p), s, r(p), 0);
  sol_err(s+1) = err(end);
  if(s > 0)
    disp('p')
    disp(p)
    disp('error')
    disp(sol_err)
    rate = (log(sol_err(1:end-1)) - log(sol_err(2:end))) / log(2);
    disp('rate')
    disp(rate)
  end
end
