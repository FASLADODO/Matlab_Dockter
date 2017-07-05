% dtemp = [-1,7];
dtemp = [-52276,0.481012658227848]

[ClassStore,RawStore] = OnlineDPPGridWeak(dtemp,ModelLoo)

%%

dist2org = dtemp - Model.origin
indexgrid = floor(dist2org ./ Model.steps) + 1
selectElement = num2cell(indexgrid)
indx=sub2ind(size(Model.DPP_Grid_Weight),selectElement{:});
Weight_temp = Model.DPP_Grid_Weight_All(selectElement{:},:)
Thresh_temp = Model.DPP_Grid_Thresh(selectElement{:},:)
Class_temp = Model.DPP_Grid_Class(selectElement{:},:,:)

Class_temp(:,:,1,:)
Class_temp(:,:,2,:)


for dd = 1:SS
    isbelow = dtemp(dd) < Thresh_temp(dd)
    dir = Class_temp(:,:,dd,:)
    classest = dir(isbelow + 1)
end