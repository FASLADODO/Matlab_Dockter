function [indx] = ComputeGridIndex(Model,Data)
%based on position and grid models compute the index in the grid

    dist2org = Data - Model.origin;
    indexgrid = floor(dist2org ./ Model.steps) + 1;
    selectElement = num2cell(indexgrid);
    indx=sub2ind(size(Model.DPP_Grid),selectElement{:});


end