function DataStruct = getClassData(Data,Labels)
%return each classes data in a seperate struct member

cslist = unique(Labels);


for ii = 1:length(cslist)
   temp = Data(Labels == cslist(ii),:);
   DataStruct{ii} = temp;
end

end