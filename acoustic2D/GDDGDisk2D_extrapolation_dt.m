clear
N = [11, 18, 25, 31, 37];
format long e
for p=1:5
  A = GDDGDisk2D_func(p, N(p), 0, 'diskhole', [], 0);
  GDDG_extrapolation_E = eig(A);
  save(sprintf('data/GDDGDisk_extrapolation_spectrum_p%d', p), 'GDDG_extrapolation_E');
end
