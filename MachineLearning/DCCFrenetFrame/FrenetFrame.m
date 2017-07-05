function [V,U,N] = FrenetFrame(Data,n)
%see Ahmidi2013 paper
%compute Frenet frame for each sequence of n points in data
%Data is a matrix with xyz columns

[NN,SS] = size(Data);

id = 1;
for tt = 1:NN-n
   %grab a window
   dtemp = Data(tt:tt+n,:);
   
   %compute the frenet frame for that window
   [vtemp,utemp,ntemp] = SingleFrenetFrame(dtemp);
   
   %stash
   V(:,id) = vtemp; %3x1
   U(:,id) = utemp; %3x1
   N(:,id) = ntemp; %3x1
   
   id = id + 1;
end