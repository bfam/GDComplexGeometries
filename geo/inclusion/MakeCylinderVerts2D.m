function MakeCylinderVerts2D(faces, ra,xo,yo)

% Function: MakeCylinderVerts2D(faces, ra, xo, yo)
% Purpose:  Make face vertices so they conform to a cylinder of radius r
%           centered at (xo,yo)
Globals2D;

NCurveFaces = size(faces,1);
for n=1:NCurveFaces

  % move vertices of faces to be curved onto circle
  k = faces(n,1); f = faces(n,2);
  v1 = EToV(k, f); v2 = EToV(k, mod(f,Nfaces)+1);

  % compute polar angles of start and end face vertices relative to circle center
  theta1 = atan2(VY(v1)-yo,VX(v1)-xo);
  theta2 = atan2(VY(v2)-yo,VX(v2)-xo);

  % move vertices onto circle
  newx1 = xo + ra*cos(theta1); newy1 = yo + ra*sin(theta1);
  newx2 = xo + ra*cos(theta2); newy2 = yo + ra*sin(theta2);

  % update mesh vertex locations
  VX(v1) = newx1; VX(v2) = newx2; VY(v1) = newy1; VY(v2) = newy2;
end

return
