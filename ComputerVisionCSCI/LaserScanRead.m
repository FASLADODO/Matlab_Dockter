% laser scan recreate

data =load('scandata.csv');


[NN,SS] = size(data);

DataAll = [];
for ff = 1:NN
   rowvals = data(ff,:); %Z height
   scanid = [1:SS]; %Y angle
   fileid = ones(1,SS)*ff; %X scan index
    
   %remove bad rows
   idr = find(rowvals == -1);
   rowvals(idr) = [];
   scanid(idr) = [];
   fileid(idr) = [];
   
   %convert
   xdata = fileid';
   ydata = scanid';
   zdata = rowvals';
   
   DataAll = [DataAll; xdata, ydata, zdata ];
end

%invert based on minimum
invertPoint = max(DataAll(:,3));
DataAll(:,3) = invertPoint - DataAll(:,3);

figure
scatter3(DataAll(:,1),DataAll(:,2),DataAll(:,3))