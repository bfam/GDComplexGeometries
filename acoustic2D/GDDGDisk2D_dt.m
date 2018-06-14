clear
N = [9, 14, 19, 24, 28];
format long e
for p = 1:5
  dt(p) = GDDGDisk2D_func(p, N(p), 0, 'diskhole');
end
disp('dt'); disp(dt)
