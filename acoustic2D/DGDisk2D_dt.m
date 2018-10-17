clear
format long e
for p = 1:5
  A = DGDisk2D_func(p, 0, 'disk');
  DG_E = eig(A);
  save(sprintf('data/DGDisk_spectrum_p%d', p), 'DG_E');
end
