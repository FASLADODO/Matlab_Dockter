function H = EntropySample(X)
%Compute entropy for a given set of categorical data
%X should be a column vector of disscrete values
%EG: X = [1,2,3,2...] or X = [1,1,1,0,0,0,1,1]

[~,SS] = size(X);

if(SS > 1)
   warning('ignoring extra columns of X');
   X = X(:,1);
end

% Get frequency table 
tab = tabulate(X);
prob = tab(:,3) / 100;
% Filter out zero-entries
prob = prob(prob~=0);
% Get entropy
H = -sum(prob .* log2(prob));


%OLD WAY 

% %get all possible categories
% cslist = unique(X);
% 
% P_X = [];
% for cc = 1:length(cslist) % loop through all categories
%    X_C = X == cslist(cc); %0 or 1 if equal to current category
%    P_X(cc) = mean(X_C); %probability of X being the current category
% end
% 
% H = Entropy(P_X);