function [vx, vy, p] = DiskExactMode2D(t, x, y, k)

p  = sqrt(2) * cos(pi * k * (x+1)/2) .* cos(pi * k * (y+1)/2) .* cos(pi * sqrt(2) * k * t/2);
vx = sin(pi * k * (x+1)/2) .* cos(pi * k * (y+1)/2) .* sin(pi * sqrt(2) * k * t/2);
vy = cos(pi * k * (x+1)/2) .* sin(pi * k * (y+1)/2) .* sin(pi * sqrt(2) * k * t/2);
