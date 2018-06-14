clear
format long e
for p = 1:5
  dt(p) = DGDisk2D_func(p, 0, 'disk');
end
disp('dt'); disp(dt)
