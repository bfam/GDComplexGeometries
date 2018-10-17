clear
N = [4, 7, 10, 14, 16];
format long e
for p=1:5
  A = GDDiskDriver2D_func(p, N(p), 0, [], 0);
  GD_extrapolation_E = eig(A);
  save(sprintf('data/GDDisk_extrapolation_spectrum_p%d', p), 'GD_extrapolation_E');
end
