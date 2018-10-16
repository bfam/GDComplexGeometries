clear
taylor_order = [4, 7, 8, 11, 12];
dtDG = zeros(5, 1);
dtGDE = zeros(5, 1);
dtGDG = zeros(5, 1);
for p = 1:5
  load(sprintf('data/box_spectrum_p%d.mat', p));
  dtDG(p) = max_Taylor_timestep(taylor_order(p), DG_E);
  dtGDE(p) = max_Taylor_timestep(taylor_order(p), GDE_E);
  dtGDG(p) = max_Taylor_timestep(taylor_order(p), GDG_E);
end

[dtDG;dtGDE;dtGDG];

% Scaled time steps
ldg = 2 / 2^2;
hdg = ldg / 2;
rdg = sqrt(2) * hdg / (1 + sqrt(2));
n = 2 * (1:5)' + 1;

Ne = [15 27 30 38 43]';
Ng = [11 24 23 30 34]';

hGDE = 2 ./ Ne;
hGDG = 2 ./ Ng;

% Data for Table 2
disp('[dtGDE./dtDG,dtGDG./dtDG]')
disp([dtGDE./dtDG,dtGDG./dtDG])

q = 2*(1:5)'+1;
NpDG_NpGDE = (4*4^2*(q+1).*(q+2) / 2) ./ (Ne+1).^2;
NpDG_NpGDG = (4*4^2* (q+1).*(q+2) / 2) ./ (Ng+1+2*(1:5)').^2;

disp('[NpDG_NpGDE, NpDG_NpGDG]')
disp([NpDG_NpGDE, NpDG_NpGDG])

disp('[NpDG_NpGDE, NpDG_NpGDG] .* [dtGDE./dtDG,dtGDG./dtDG]')
disp([NpDG_NpGDE, NpDG_NpGDG] .* [dtGDE./dtDG,dtGDG./dtDG])

% Data for Table 3
disp('[n .* dtDG / (2*rdg), dtGDE ./ hGDE, dtGDG ./ hGDG]')
disp([n .* dtDG / (2*rdg), dtGDE ./ hGDE, dtGDG ./ hGDG])

