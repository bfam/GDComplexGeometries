clear

format short e

Ne = [15 27 30 38 43];
Ng = [11 24 23 30 34];
for p = 1:5
  DG_A = DGbox_func(p);
  DG_E = eig(DG_A);

  GDE_A = GDbox_func(p, 0, Ne(p));
  GDE_E = eig(GDE_A);

  GDG_A = GDbox_func(p, p, Ng(p));
  GDG_E = eig(GDG_A);
  save(sprintf('data/box_spectrum_p%d', p), 'DG_E',  'GDE_E',  'GDG_E')
end
