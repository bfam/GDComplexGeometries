function [mapB, mapM] = BuildBCMaps2D_cubature(intC, BCType, vmapM)

bct    = BCType';
bnodes = ones(intC, 1)*bct(:)';
bnodes = bnodes(:);

mapB = find(bnodes == 1);
mapM = find(bnodes > 1);
