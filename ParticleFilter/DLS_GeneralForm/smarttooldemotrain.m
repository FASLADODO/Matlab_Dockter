%smarttooldemo train

clear all

%Area of the load cell we think maybe
AREAZ = 1.86325e-04;

hard_id = 1;
soft_id = 2;

fileID1 = fopen('softdata.csv');
Dtemp1 = textscan(fileID1,'%f%f%f%f%f%f', 'delimiter',',');
Dall1 = double(cell2mat(Dtemp1));

FullData{hard_id}.t = Dall1(:,1);
FullData{hard_id}.ind = Dall1(:,2);
FullData{hard_id}.stress = Dall1(:,3);
FullData{hard_id}.stressdot = Dall1(:,4);
FullData{hard_id}.stressdotdot = Dall1(:,5);
FullData{hard_id}.angle = Dall1(:,6);
FullData{hard_id}.angledot = Dall1(:,7);
FullData{hard_id}.angledotdot = Dall1(:,8);
FullData{hard_id}.graspenum = Dall1(:,9);


%pre allocate
for ii = 1:max(FullData{hard_id}.ind)
    Data_Hard{ii}.t = [];
    Data_Hard{ii}.ind = [];
    Data_Hard{ii}.stress = [];
    Data_Hard{ii}.stressdot = [];
    Data_Hard{ii}.stressdotdot = [];
    Data_Hard{ii}.angle = [];
    Data_Hard{ii}.angledot = [];
    Data_Hard{ii}.angledotdot = [];
end
%fill in each grasp
for jj = 1:max(FullData{hard_id}.ind)
    if(FullData{hard_id}.graspenum(jj) == 2)
        Data_Hard{FullData{hard_id}.ind(jj)}.t = [Data_Hard{FullData{hard_id}.ind(jj)}.t; FullData{hard_id}.t(jj)];
        Data_Hard{FullData{hard_id}.ind(jj)}.stress = [Data_Hard{FullData{hard_id}.ind(jj)}.stress; FullData{hard_id}.stress(jj)];
        Data_Hard{FullData{hard_id}.ind(jj)}.stressdot = [Data_Hard{FullData{hard_id}.ind(jj)}.stressdot; FullData{hard_id}.stressdot(jj)];
        Data_Hard{FullData{hard_id}.ind(jj)}.stressdotdot = [Data_Hard{FullData{hard_id}.ind(jj)}.stressdotdot; FullData{hard_id}.stressdotdot(jj)];
        Data_Hard{FullData{hard_id}.ind(jj)}.angle = [Data_Hard{FullData{hard_id}.ind(jj)}.angle; FullData{hard_id}.angle(jj)];
        Data_Hard{FullData{hard_id}.ind(jj)}.angledot = [Data_Hard{FullData{hard_id}.ind(jj)}.angledot; FullData{hard_id}.angledot(jj)];
        Data_Hard{FullData{hard_id}.ind(jj)}.angledotdot = [Data_Hard{FullData{hard_id}.ind(jj)}.angledotdot; FullData{hard_id}.angledotdot(jj)];
        
    end
end

fileID2 = fopen('softdata.csv');
Dtemp2 = textscan(fileID2,'%f%f%f%f%f%f', 'delimiter',',');
Dall2 = double(cell2mat(Dtemp2));

FullData{soft_id}.t = Dall2(:,1);
FullData{soft_id}.ind = Dall2(:,2);
FullData{soft_id}.stress = Dall2(:,3);
FullData{soft_id}.stressdot = Dall2(:,4);
FullData{soft_id}.stressdotdot = Dall2(:,5);
FullData{soft_id}.angle = Dall2(:,6);
FullData{soft_id}.angledot = Dall2(:,7);
FullData{soft_id}.angledotdot = Dall2(:,8);

%pre allocate
for ii = 1:max(FullData{soft_id}.ind)
    Data_Soft{ii}.t = [];
    Data_Soft{ii}.ind = [];
    Data_Soft{ii}.stress = [];
    Data_Soft{ii}.stressdot = [];
    Data_Soft{ii}.stressdotdot = [];
    Data_Soft{ii}.angle = [];
    Data_Soft{ii}.angledot = [];
    Data_Soft{ii}.angledotdot = [];
end
%fill in each grasp
for jj = 1:max(FullData{soft_id}.ind)
    if(FullData{soft_id}.graspenum(jj) == 2)
        Data_Soft{FullData{soft_id}.ind(jj)}.t = [Data_Soft{FullData{soft_id}.ind(jj)}.t; FullData{soft_id}.t(jj)];
        Data_Soft{FullData{soft_id}.ind(jj)}.stress = [Data_Soft{FullData{soft_id}.ind(jj)}.stress; FullData{soft_id}.stress(jj)];
        Data_Soft{FullData{soft_id}.ind(jj)}.stressdot = [Data_Soft{FullData{soft_id}.ind(jj)}.stressdot; FullData{soft_id}.stressdot(jj)];
        Data_Soft{FullData{soft_id}.ind(jj)}.stressdotdot = [Data_Soft{FullData{soft_id}.ind(jj)}.stressdotdot; FullData{soft_id}.stressdotdot(jj)];
        Data_Soft{FullData{soft_id}.ind(jj)}.angle = [Data_Soft{FullData{soft_id}.ind(jj)}.angle; FullData{soft_id}.angle(jj)];
        Data_Soft{FullData{soft_id}.ind(jj)}.angledot = [Data_Soft{FullData{soft_id}.ind(jj)}.angledot; FullData{soft_id}.angledot(jj)];
        Data_Soft{FullData{soft_id}.ind(jj)}.angledotdot = [Data_Soft{FullData{soft_id}.ind(jj)}.angledotdot; FullData{soft_id}.angledotdot(jj)];
        
    end
end

%count number of grasps
num_Hard_Grasps = length(Data_Hard);
num_Soft_Grasps = length(Data_Soft);


data_all = [];
data_all{1}.state = [];
data_all{1}.input = [];
data_all{2}.state = [];
data_all{2}.input = [];

%Get All Params for TLS
for ii = 1:num_Hard_Grasps
    
    data_all{1}.state = [data_all{1}.state; Data_Hard(ii).strainDotDot.*(Data_Hard(ii).firstTouch.^2), Data_Hard(ii).strainDot.*Data_Hard(ii).firstTouch, Data_Hard(ii).strain, Data_Hard(ii).strain.^2, Data_Hard(ii).strain.^3];
    data_all{1}.input = [data_all{1}.input; Data_Hard(ii).stress];

end
Params_Hard_All = pinv(data_all{1}.state)*data_all{1}.input;


Residual_Hard=data_all{1}.input-data_all{1}.state*Params_Hard_All;


for ii = 1:num_Soft_Grasps
    data_all{2}.state = [data_all{2}.state; Data_Soft(ii).strainDotDot.*(Data_Soft(ii).firstTouch.^2), Data_Soft(ii).strainDot.*Data_Soft(ii).firstTouch, Data_Soft(ii).strain, Data_Soft(ii).strain.^2, Data_Soft(ii).strain.^3];
    data_all{2}.input  = [data_all{2}.input; Data_Soft(ii).stress];
end


Params_Soft_All = pinv(data_all{2}.state)*data_all{2}.input;

Residual_Soft=data_all{2}.input-data_all{2}.state*Params_Soft_All;