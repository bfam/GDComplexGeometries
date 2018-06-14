function [fout] = RefineUniformField2D(RefineLevel, f)

Globals2D;

[rout, sout] = RefineUniformRS2D(r,s);
IM =  InterpMatrix2D(rout, sout);

fout = f;
for level = 1:RefineLevel
  fout = IM*fout;
  fout = reshape(fout, Np, 4*size(fout,2));
end

end
