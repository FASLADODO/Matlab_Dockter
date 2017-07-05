function [Data,Labels] = CreateData(means,sigmas,nn)
%Creates random data

classes = length(means);
dims = length(means{1});
sdims = length(sigmas{1});

if(dims ~= sdims)
    error('size of means needs to equal size of sigmas')
end

if(length(nn) == 1)
   nn = ones(classes,1)*nn(1); 
end

Data = [];
Labels = [];

for ii = 1:classes
    dt = mvnrnd(means{ii},sigmas{ii},nn(ii));
    Data = [Data; dt];
    Labels = [Labels; ones(nn(ii),1)*ii ];
end

end