function MakeInclusionsVerts2D()

Globals2D;

cyl.bc = [5, 6, 7, 8];
cyl.r  = [2, 2.5, 1.5, 1];
cyl.x  = [2, -2, -3, 2.5];
cyl.y  = [2, -2, 2.5, -3];
for c = 1:length(cyl.bc)
  [k, f] = find(BCType == cyl.bc(c));
  if(~isempty(k))
    outfaces = [k,f];
    MakeCylinderVerts2D(outfaces, cyl.r(c), cyl.x(c), cyl.y(c));
  end
end

end
