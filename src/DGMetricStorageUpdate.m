function C = DGMetricStorageUpdate(C)
Globals2D;

C.W = C.W ./ C.J;
C.J = C.V * J;

C.rx = (C.V * (J .* rx)) ./ C.J;
C.ry = (C.V * (J .* ry)) ./ C.J;
C.sx = (C.V * (J .* sx)) ./ C.J;
C.sy = (C.V * (J .* sy)) ./ C.J;
C.WJI = C.W ./ C.J;
C.W = C.W .* C.J;
