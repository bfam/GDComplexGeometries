clear
N = [4, 6, 7, 11, 13];
format long e
for p=1:5
  dt(p) = GDDiskDriver2D_func(p, N(p), 0);
end
disp('dt'); disp(dt)
