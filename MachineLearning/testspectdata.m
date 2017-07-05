%test spect data with naive bayes

fileID = fopen('spectdata.txt');
Dtemp = textscan(fileID,'%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d', 'delimiter',',');
%titles = {DIAGNOSIS,F1,F2,F3,F4,F5,F6,F7,F8,F9,F10,F11,F12,F13,F14,F15,F16,F17,F18,F19,F20,F21,F22};


Dall = double(cell2mat(Dtemp));

Data = Dall(:,2:end);
Y = Dall(:,1);

[NN,SS] = size(Data);

%% get params for NB

%get all zeros in Y
notY = 1 - Y;

V = 2; %laplace smoothing parameter

phi_y_1 = sum(Y)/NN; %p(y=1)
phi_y_0 = 1-phi_y_1; %p(y=0)
phi_1 = []; %p(x=1|y=1)
phi_0 = []; %p(x=1|y=0)
%loop through and get parameters
for jj = 1:SS
   
    for kk = 0:1
        phi_1(kk+1,jj) = (sum( double(Data(:,jj)==kk) & Y) + 1) / (sum(Y) + V);
        phi_0(kk+1,jj) = (sum( double(Data(:,jj)==kk) & notY) + 1) / (sum(notY)+ V);
    end
end

%remove any zeros or ones cuz thats bogus
% for jj = 1:SS
%     for kk = 0:1
%         if(phi_1(kk+1,jj) == 0)
%             phi_1(kk+1,jj) = phi_1(kk+1,jj) + 0.001;
%         elseif(phi_1(kk+1,jj) == 1)
%             phi_1(kk+1,jj) = phi_1(kk+1,jj) - 0.001;
%         end
%         if(phi_0(kk+1,jj) == 0)
%             phi_0(kk+1,jj) = phi_0(kk+1,jj) + 0.001;
%         elseif(phi_0(kk+1,jj) == 1)
%             phi_0(kk+1,jj) = phi_0(kk+1,jj) - 0.001;
%         end
%     end
% end



sum(phi_0)
sum(phi_1)


%% Test classify for NB


Probz = [];
for ii = 1:NN
    prod_1 = 1;
    prod_0 = 1;
    for jj = 1:SS
        prod_1 = prod_1*phi_1(Data(ii,jj)+1,jj);
        prod_0 = prod_0*phi_0(Data(ii,jj)+1,jj);
    end
    
    p_1_x = (prod_1*phi_y_1) /( prod_1*phi_y_1 + prod_0*phi_y_0);
    p_0_x = (prod_0*phi_y_0) /( prod_1*phi_y_1 + prod_0*phi_y_0);
    
    [mx, estclass] = max([p_0_x,p_1_x]);
    Probz = [Probz; Y(ii), estclass-1, p_0_x, p_1_x];
end


%This thing is crazy good
classify = Probz(:,1) == Probz(:,2);
accuracy = sum(classify)/length(classify)