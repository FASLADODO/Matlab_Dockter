%test naive

fileID = fopen('nurserydata.txt');
D=textscan(fileID,'%s%s%s%s%s%s%s%s%s', 'delimiter',',');

nn = length(D{1});


%possible string values for each topic
parents = {'usual', 'pretentious', 'great_pret'}; %1
has_nurs = {'proper', 'less_proper', 'improper', 'critical', 'very_crit'}; %2
form = {'complete', 'completed', 'incomplete', 'foster'}; %3
children = {'1', '2', '3', 'more'}; %4
housing = {'convenient', 'less_conv', 'critical'}; %5
finance= {'convenient', 'inconv'}; %6
social = {'nonprob', 'slightly_prob', 'problematic'}; %7
health = {'recommended', 'priority', 'not_recom'}; %8
result = {'not_recom', 'recommend', 'very_recom', 'priority', 'spec_prior'}; %9

%get it in binary
Data = [];
Y = [];
for ii = 1:nn
    tempd = [];
    
    tempy = double(strcmp(D{9}{ii},result));
    
    if(tempy(4)) %'not_recom'
        tempd = [tempd, find( double(strcmp(D{1}{ii},parents)) == 1)];
        tempd = [tempd, find( double(strcmp(D{2}{ii},has_nurs)) == 1)];
        tempd = [tempd, find( double(strcmp(D{3}{ii},form)) == 1)];
        tempd = [tempd, find( double(strcmp(D{4}{ii},children)) == 1)];
        tempd = [tempd, find( double(strcmp(D{5}{ii},housing)) == 1)];
        tempd = [tempd, find( double(strcmp(D{6}{ii},finance)) == 1)];
        tempd = [tempd, find( double(strcmp(D{7}{ii},social)) == 1)];
        %tempd = [tempd, find( double(strcmp(D{8}{ii},health)) == 1)];
        Data = [Data; tempd];
        Y = [Y; 0];
    elseif( tempy(5) ) %'priority', 'spec_prior'
        tempd = [tempd, find( double(strcmp(D{1}{ii},parents)) == 1)];
        tempd = [tempd, find( double(strcmp(D{2}{ii},has_nurs)) == 1)];
        tempd = [tempd, find( double(strcmp(D{3}{ii},form)) == 1)];
        tempd = [tempd, find( double(strcmp(D{4}{ii},children)) == 1)];
        tempd = [tempd, find( double(strcmp(D{5}{ii},housing)) == 1)];
        tempd = [tempd, find( double(strcmp(D{6}{ii},finance)) == 1)];
        tempd = [tempd, find( double(strcmp(D{7}{ii},social)) == 1)];
        %tempd = [tempd, find( double(strcmp(D{8}{ii},health)) == 1)];
        Data = [Data; tempd];
        Y = [Y; 1];
    end
end

Variates = max(Data);

%% Get parameters

[NN,SS] = size(Data);

%get all zeros in Y
notY = 1 - Y;


phi_y_1 = sum(Y)/NN; %p(y=1)
phi_y_0 = 1-phi_y_1; %p(y=0)
phi_1 = []; %p(x=1|y=1)
phi_0 = []; %p(x=1|y=0)
%loop through and get parameters
for jj = 1:SS
    Variates(jj)
    for kk = 1:Variates(jj)
        phi_1(kk,jj) = sum( double(Data(:,jj)==kk) & Y) / sum(Y);
        phi_0(kk,jj) = sum( double(Data(:,jj)==kk) & notY) / sum(notY);
    end
end

%remove any zeros or ones cuz thats bogus
% for jj = 1:SS
%     if(phi_1_1(jj) == 0)
%         phi_1_1(jj) = phi_1_1(jj) + 0.001;
%     elseif(phi_1_1(jj) == 1)
%         phi_1_1(jj) = phi_1_1(jj) - 0.001;
%     end
%     if(phi_1_0(jj) == 0)
%         phi_1_0(jj) = phi_1_0(jj) + 0.001;
%     elseif(phi_1_0(jj) == 1)
%         phi_1_0(jj) = phi_1_0(jj) - 0.001;
%     end
% end

sum(phi_0)
sum(phi_1)

%% Test classify


Probz = [];
for ii = 1:NN
    prod_1 = 1;
    prod_0 = 1;
    for jj = 1:SS
        prod_1 = prod_1*phi_1(Data(ii,jj),jj);
        prod_0 = prod_0*phi_0(Data(ii,jj),jj);
    end
    
    p_1_x = (prod_1*phi_y_1) /( prod_1*phi_y_1 + prod_0*phi_y_0);
    p_0_x = (prod_0*phi_y_0) /( prod_1*phi_y_1 + prod_0*phi_y_0);
    
    [mx, estclass] = max([p_0_x,p_1_x]);
    Probz = [Probz; Y(ii), estclass-1, p_0_x, p_1_x];
end


%This thing is crazy good
classify = Probz(:,1) == Probz(:,2);
accuracy = sum(classify)/length(classify)




%%

K1 = Data(:,1) == 1;
K2 = Data(:,1) == 2;
K3 = Data(:,1) == 3;

rat1 = sum(K1)/NN
rat2 = sum(K2)/NN
rat3 = sum(K3)/NN




