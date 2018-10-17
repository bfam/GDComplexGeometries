clear
N = [4, 6, 7, 11, 13];
format long e
for p=1:5
  A = GDDiskDriver2D_func(p, N(p), 0);
  GD_E = eig(A);
  save(sprintf('data/GDDisk_spectrum_p%d', p), 'GD_E');
end
