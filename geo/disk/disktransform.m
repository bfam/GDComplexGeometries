% setup the curved blocks for a disk
function z = disk_transform(r1, r2, block, var)

d = 1/2;

if block == 1
  x = r1 + 0*r2;
  y = r2 + 0*r1;
elseif block == 2
  a = (1+r1)/2;
  R = 3*(1-a) + a;
  Q = pi - r2 * pi /4;

  x0 = -((3*sqrt(2)/2)*(1-a) + a);
  y0 = -x0 .* r2;

  xt = R .* cos(Q);
  yt = R .* sin(Q);

  x = xt .* (1-a) + x0 .* a;
  y = yt .* (1-a) + y0 .* a;

elseif block == 3
  a = (1+r1)/2;
  R = 3*a + (1-a);
  Q = r2 * pi /4;

  x0 = (3*sqrt(2)/2)*a + (1-a);
  y0 = x0 .* r2;

  yt = R .* sin(Q);
  xt = R .* cos(Q);

  x = xt .* a + x0 .* (1-a);
  y = yt .* a + y0 .* (1-a);

elseif block == 4
  b = (1+r2)/2;
  R = 3*(1-b) + b;
  Q = 3*pi/2 + r1 * pi /4;

  y0 = -((3*sqrt(2)/2)*(1-b) + b);
  x0 = -y0 .* r1;

  xt = R .* cos(Q);
  yt = R .* sin(Q);

  x = xt .* (1-b) + x0 .* b;
  y = yt .* (1-b) + y0 .* b;

elseif block == 5
  b = (r2+1)/2;
  R = 3*b + (1-b);
  Q = r1 * pi /4;

  y0 = (3*sqrt(2)/2)*b + (1-b);
  x0 = y0 .* r1;

  yt = R .* cos(Q);
  xt = R .* sin(Q);

  y = yt .* b + y0 .* (1-b);
  x = xt .* b + x0 .* (1-b);
end

if var == 1
  z = x/3;
else
  z = y/3;
end
