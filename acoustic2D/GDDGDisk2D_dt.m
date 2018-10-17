clear
N = [9, 14, 19, 24, 28];
format long e
for p = 1:5
  A = GDDGDisk2D_func(p, N(p), 0, 'diskhole');
  GDDG_E = eig(A);
  save(sprintf('data/GDDGDisk_spectrum_p%d', p), 'GDDG_E');
end
