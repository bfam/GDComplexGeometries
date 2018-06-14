function [rur, rus] = RefineUniformRS2D(r, s)

% Function: RefineUniformRS2D(r, s)
% Purpose:  Provide the r and s coordiants for one level of refinement

rur = [ r/2 - 1/2;
        r/2 + 1/2;
        r/2 - 1/2;
       -s/2 - 1/2 ];

rus = [ s/2 - 1/2;
        s/2 - 1/2;
        s/2 + 1/2;
        r/2 + s/2 ];

end
