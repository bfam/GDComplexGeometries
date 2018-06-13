clear
N = [4, 7, 10, 14, 16];
format long e
for p=1:5
  dt(p) = GDDiskDriver2D_func(p, N(p), 0, [], 0);
end
disp('dt'); disp(dt)
