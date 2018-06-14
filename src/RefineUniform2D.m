function RefineUniform2D(RefineLevel, MoveVerts)

% Function: RefineUniform2D(RefineLevel, MoveVerts)
% Purpose:  Uniformly refines a mesh RefineLevel times and calls MoveVerts each
%           iteration.


Globals2D;

Nv = length(VX);

for level=1:RefineLevel
  newK = 4*K;
  newEToV = zeros(newK, 3);
  newBCType = zeros(newK, 3);

  newVX = zeros(1,3*newK);
  newVY = zeros(1,3*newK);
  newVX(1:Nv) = VX;
  newVY(1:Nv) = VY;

  % Refine mesh
  Vmap = sparse(Nv,Nv);
  v = Nv;
  for k = 1:K
    v1 = EToV(k, 1); v2 = EToV(k, 2); v3 = EToV(k, 3);

    v4 = Vmap(v1, v2); v5 = Vmap(v2, v3); v6 = Vmap(v1, v3);
    if v4 == 0
      v = v + 1;
      v4 = v;
      Vmap(v1, v2) = v4;
      Vmap(v2, v1) = v4;
      newVX(v) = (VX(v1) + VX(v2))/2;
      newVY(v) = (VY(v1) + VY(v2))/2;
    end
    if v5 == 0
      v = v + 1;
      v5 = v;
      Vmap(v2, v3) = v5;
      Vmap(v3, v2) = v5;
      newVX(v) = (VX(v2) + VX(v3))/2;
      newVY(v) = (VY(v2) + VY(v3))/2;
    end
    if v6 == 0
      v = v + 1;
      v6 = v;
      Vmap(v1, v3) = v6;
      Vmap(v3, v1) = v6;
      newVX(v) = (VX(v1) + VX(v3))/2;
      newVY(v) = (VY(v1) + VY(v3))/2;
    end

    newEToV(4*(k-1)+1, 1) = v1; newEToV(4*(k-1)+1, 2) = v4; newEToV(4*(k-1)+1, 3) = v6;
    newEToV(4*(k-1)+2, 1) = v4; newEToV(4*(k-1)+2, 2) = v2; newEToV(4*(k-1)+2, 3) = v5;
    newEToV(4*(k-1)+3, 1) = v6; newEToV(4*(k-1)+3, 2) = v5; newEToV(4*(k-1)+3, 3) = v3;
    newEToV(4*(k-1)+4, 1) = v4; newEToV(4*(k-1)+4, 2) = v5; newEToV(4*(k-1)+4, 3) = v6;

    newBCType(4*(k-1)+1, 1) = BCType(k, 1); newBCType(4*(k-1)+1, 2) = 0;            newBCType(4*(k-1)+1, 3) = BCType(k, 3);
    newBCType(4*(k-1)+2, 1) = BCType(k, 1); newBCType(4*(k-1)+2, 2) = BCType(k, 2); newBCType(4*(k-1)+2, 3) = 0;
    newBCType(4*(k-1)+3, 1) = 0;            newBCType(4*(k-1)+3, 2) = BCType(k, 2); newBCType(4*(k-1)+3, 3) = BCType(k, 3);
    newBCType(4*(k-1)+4, 1) = 0;            newBCType(4*(k-1)+4, 2) = 0;            newBCType(4*(k-1)+4, 3) = 0;
  end
  newNv = v;
  newVX = newVX(1:newNv);
  newVY = newVY(1:newNv);

  Nv = newNv;
  VX = newVX;
  VY = newVY;
  K = newK;
  EToV = newEToV;
  BCType = newBCType;

  MoveVerts();
end

end
