function [W, C] = LDASimple(X,Y,Priors)
%super basic LDA code
% Use:
% [W,C] = LDA(X,Y,Priors)
%
% W     = LDA coefficients, each row is for a class
% C     = constants, each row is for a class
% X   = predictor data (variables in columns, observations in rows)
% Y  = target variable (class labels) :1,2,3...
% Priors  = vector of prior probabilities (optional)

%Sizes
[NN,SS] = size(X);

%unique class values
ClassLabel = unique(Y);
g = length(ClassLabel);

%placeholder
mu = [];

%pooled covariance
PCov=zeros(SS,SS);

%loop through each unique class
for ii = 1:g 
    %get all samples for current class
   idx = find(Y == ClassLabel(ii));
   groupX =  X(idx,:);
   [ni,si] = size(groupX);
   
   %get means
   mu(ii,:) = mean(groupX,1);
   
   %add to pooled covariance (cov() takes care of shifting)
   PCov = PCov + ( (ni-1)/(NN-g)) .* cov( groupX ); 
   
   %probabilities
   PriorProb(ii) = ni/NN;
end

%Use actual priors for probabilities if provided
if(nargin == 3)
    PriorProb = Priors;
end


% Loop over classes to calculate linear discriminant coefficients
for ii = 1:g
    % Intermediate calculation for efficiency
    % This replaces:  mu(i,:) * inv(PCov)
    Temp = mu(ii,:) / PCov;
    
    % Constant
    C(ii,:) = -0.5 * Temp * mu(ii,:)' + log(PriorProb(ii));
    
    % Linear coefficients
    W(ii,:) = Temp;
end

Map = ClassLabel; %useful with MapValues

%returns W and C

end