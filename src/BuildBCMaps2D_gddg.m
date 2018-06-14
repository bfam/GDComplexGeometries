% build the maps for the boundary conditions and/or mortar for the nodal-dg
% codes
function BuildBCMaps2D_gddg(bc)

if nargin < 1
  bc = 1;
end


Globals2D_gddg;
Globals2D;
mor = unique(BCType);
mor = setdiff(mor, 0);
mor = setdiff(mor, bc);

bct    = BCType';
bnodes = ones(Nfp, 1)*bct(:)';
bnodes = bnodes(:);

mapB = [];
for k = 1:length(bc)
  mapB = find(bnodes == bc(k));
end
mapB = sort(mapB);
vmapB = vmapM(mapB);

mapMor = {};
vmapMor = {};

for k = 1:length(mor)
  mapMor{k} = find(bnodes == mor(k));
  mapMor{k} = reshape(mapMor{k}, Nfp, length(mapMor{k})/Nfp);
  vmapMor{k} = vmapM(mapMor{k});
end
