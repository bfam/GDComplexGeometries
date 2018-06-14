clear
N = [11, 18, 25, 31, 37];
format long e
for p=1:5
  dt(p) = GDDGDisk2D_func(p, N(p), 0, 'diskhole', [], 0);
end
disp('dt'); disp(dt)
