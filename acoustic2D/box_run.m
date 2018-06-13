clear

format short e

k  = [1,  4, 7, 12, 15];
Ne = [15 27 30 38 43];
Ng = [11 24 23 30 34];

dtDG  = [6.4396e-03   4.3856e-03   3.3445e-03   2.7091e-03   2.2770e-03];
dtGDG = [1.6360e-02   1.0521e-02   7.2008e-03   5.4346e-03   4.5745e-03];
dtGDE = [1.1878e-02   8.0147e-03   4.8555e-03   3.7588e-03   3.6733e-03];

% 1
% 4
% 3
% 4
% 5

for p = 1:5
  out = DGbox_func(p, dtDG(p)); dtDG(p) = out(1); dt0DG(p) = out(2);
  out = DGbox_func(p,  k(p)  ,  dtDG(p) / 10); errDG(p) = out(end);
  out = DGbox_func(p,  k(p)+1,  dtDG(p) / 10); errDG1(p) = out(end);

  out = GDbox_func(p, 0, Ne(p), dtGDE(p)); dtGDE(p) = out(1); dt0GDE(p) = out(2);
  out = GDbox_func(p, 0, Ne(p)  , k(p), dtGDE(p) / 10); errGDE(p) = out(end);
  out = GDbox_func(p, 0, Ne(p)-1, k(p), dtGDE(p) / 10); errGDE1(p) = out(end);

  out = GDbox_func(p, p, Ng(p), dtGDG(p)); dtGDG(p) = out(1); dt0GDG(p) = out(2);
  out = GDbox_func(p, p, Ng(p)  , k(p), dtGDG(p) / 10); errGDG(p) = out(end);
  out = GDbox_func(p, p, Ng(p)-1, k(p), dtGDG(p) / 10); errGDG1(p) = out(end);
end

disp('errDG < 1e-3')
disp(errDG(1:p) < 1e-3)
disp('errDG1 > 1e-3')
disp(errDG1(1:p) > 1e-3)

disp('errDG > errGDE')
disp(errDG(1:p) > errGDE)
disp('errDG < errGDE1')
disp(errDG(1:p) < errGDE1)

disp('errDG > errGDG')
disp(errDG(1:p) > errGDG)
disp('errDG < errGDG1')
disp(errDG(1:p) < errGDG1)

disp('[errDG; errDG1]')
disp([errDG; errDG1])
disp('[errGDE; errGDE1]')
disp([errGDE; errGDE1])
disp('[errGDG; errGDG1]')
disp([errGDG; errGDG1])

disp('[dtGDE./dtDG;dtGDG./dtDG]')
disp([dtGDE./dtDG;dtGDG./dtDG])
disp('[dtDG;dtGDE;dtGDG]')
disp([dtDG;dtGDE;dtGDG])
disp('[dt0DG;dt0GDE;dt0GDG]')
disp([dt0DG;dt0GDE;dt0GDG])
disp('[dtDG;dtGDE;dtGDG]./[dt0DG;dt0GDE;dt0GDG]')
disp([dtDG;dtGDE;dtGDG]./[dt0DG;dt0GDE;dt0GDG])

q = 2*(1:5)+1;
disp('(4*4^2* (q+1).*(q+2) / 2) ./ (Ne+1).^2')
disp((4*4^2* (q+1).*(q+2) / 2) ./ (Ne+1).^2)
disp('(4*4^2* (q+1).*(q+2) / 2) ./ (Ng+1+2*(1:5)).^2')
disp((4*4^2* (q+1).*(q+2) / 2) ./ (Ng+1+2*(1:5)).^2)

%{
ldg = 2 / 2^2;
hdg = ldg / 2;
rdg = sqrt(2) * hdg / (1 + sqrt(2));
n = 2 * (1:5) + 1;

hGDE = 2 ./ Ne;
hGDG = 2 ./ Ng;

disp([n .* dtDG / rdg; dtGDE ./ hGDE; dtGDG ./ hGDG])
%}
