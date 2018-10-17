% script to run curved box test problem with extrapolation method
% NOTE: p should be defined before running this script
format long e
sol_err = [];
for s = 0:4
  [T, err] = GDbox_curved_func(p, 2^s * [15 20 20 15], 15, 0);
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

disp('FINAL: GDBox_curved_extrapolation')
disp('p')
disp(p)
disp('error')
disp(sol_err)
rate = (log(sol_err(1:end-1)) - log(sol_err(2:end))) / log(2);
disp('rate')
disp(rate)
