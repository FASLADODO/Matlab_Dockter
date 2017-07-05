function [Thresh,bestscore,IG] = DecisionStumpBasic(X,W,Labels,Direction,XSort,DX,stepratio)
%Computes optimal threshold for X numerical data to classify with
%X is a single column of data
%W are the current weights for each samples
%Labels are class labels for each sample
%Optional:
%Direction can specify the direction that threshold classifies 0/1
%XSort if X has already been sorted
%DX if step size is already known
%stepratio: default is 0.5 ie X + DX*0.5
%Returns:
%Thresh: optimal threshold
%best score is the sum of weights of incorrectly classified data
%IG: information gain from classify with this best threshold
%http://www.autonlab.org/_media/tutorials/dtree18.pdf

[NN,SS] = size(X);
if(SS > 1)
   warning('X has too many columns, using only the first')
   X = X(:,1);
end

%Fill in missing info if minimal function is called
if(nargin < 7)
    if(nargin == 6)
        stepratio = 0.5;
    end
    cslist = unique(Labels);

    if(nargin < 4)
        %figure out which class has lower mean
        for ii = 1:length(cslist)
            dtemp = X(Labels == cslist(ii),:);
            mu(ii) = mean(dtemp);
            sigma(ii) = std(dtemp);
        end
        %figure out which class has the lowest mean
        [~,muid] = sort(mu);
        Direction = cslist(muid);
        %figure out which mean has which standard deviation
        sigmu = sigma(muid);
        %change the step ratio based on std
        grad = ((min(sigma)/max(sigma)) - 1);
        stepratio = 0.5 + sign(sigmu(2)-sigmu(1))*grad*0
        %grad = (SIG(2) / SIG(1)) .* sign(sigmu(1) - sigmu(2));
        %stepratio = (atan(0.5*grad) / pi) + 0.5
    end

    %sort em
    XSort = sort(X);

    %steps
    DX = diff(XSort);
end


%initialize it
Thresh = 0;
bestscore = Inf;
IG = [];
IG_Temp = 0;

%loop thorugh all possible thresholds
for ii = 1:NN-1
    %try a threshold
    threshtemp = XSort(ii,:) + (DX(ii,:)*stepratio);
    isabove = X >= threshtemp;
    
    %see whats incorrect above and below threshold
    wrong = Direction(isabove+1) ~= Labels;
    
    %see which class that is
    ClassEst = Direction(isabove+1);
    
    %compute information gain
    %IG_Temp = InformationGain(Labels,ClassEst);
    
    %compute weights of missclassified
    score = sum(W(wrong));
    
    %should we keep it
    if(score < bestscore)
       Thresh = threshtemp;
       bestscore = score;
       IG = IG_Temp;
    end
end


end