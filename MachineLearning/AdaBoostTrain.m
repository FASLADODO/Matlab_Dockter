function [Tree] = AdaBoostTrain(X,Labels,Rounds,ploton)
%Train an Adaboost classifier from scratch
%Using the Mohri Book derivation
%Rod Dockter

%Data = matrix for training (Rows as Samples, Columns as Dimensions)
%Labels = {-1,1} class membership for each point in Data
%Rounds = 3 number of boosting iterations
%ploton = 1 to see the weak learners

%get data sizes
[NN,SS] = size(X);
[Ni,~] = size(Labels);
if(Ni ~= NN)
   error('number of data points must match number of labels') 
end

%get the class labels
cslist = unique(Labels);
if( sum(ismember(cslist,[-1,1])) ~= 2 )
    error('labels must be {-1,1}') 
end

%defaults
if(nargin < 4)
   ploton = false; 
end

%weights
Weights = [ones(NN,1)/NN];

%sort in each dimension
for jj = 1:SS
    for ii = 1:length(cslist)
        dtemp = X(Labels == cslist(ii),jj);
        mu(ii) = mean(dtemp);
    end
    %figure out which class has the lowest mean
    [~,muid] = sort(mu);
    Direction{jj} = cslist(muid);
    %sort em
    XSort{jj} = sort(X(:,jj));
    DX{jj} = diff(XSort{jj});
end

%for plots ignore
scale = 1000;

%This is the struct where everything goes
Tree = [];
Tree.Rounds = Rounds;
Tree.WeightInit = Weights;

if(ploton && SS >= 2)
    fig=figure; 
end

% for adaboost
for idt = 1:Tree.Rounds
    fprintf('round %i of %i',idt,Tree.Rounds);
    
    %initialize test variables
    score = Inf;
    threshSave = 0;
    igsave = 0;
    dirsave = [];
    dimsave = 1;
    
    %test each dimension
    for jj = 1:SS
        [Thresh,bestscore,IG] = DecisionStumpBasic(X(:,jj),Weights,Labels,Direction{jj},XSort{jj},DX{jj});
            
        %check if this parrallel plane results in the lowest error yet
        if (bestscore < score)
           score = bestscore;
           threshSave = Thresh;
           dirsave = Direction{jj};
           dimsave = jj;
           igsave = IG;
        end
    end
    %save the best decision tree and corresponding info
    Tree.Layer{idt}.Threshold = threshSave;
    Tree.Layer{idt}.Direction = dirsave;
    Tree.Layer{idt}.Dimension = dimsave;
    Tree.Layer{idt}.Accuracy = score;
    Tree.Layer{idt}.InformationGain = igsave;
    
    %get class prediction
    ClassEstTemp = AdaBoostDecision(X,Tree,idt);

    %update weights based on this learner
    [Weights,alpha] = AdaBoostUpdate(Weights,Labels,ClassEstTemp);
    Tree.Layer{idt}.Alpha = alpha; %save alpha for this weak learner
    Tree.Layer{idt}.Weight = Weights; %save resultant weights

    %this will be annoying for big rounds
    if(ploton && SS >= 2)
        clf(fig)
        hax=axes; 
        otherdim = SS + 1 - dimsave;
        d1t = X(ClassEstTemp == dirsave(1),:);
        d2t = X(ClassEstTemp == dirsave(2),:);
        w1t = Weights(ClassEstTemp == dirsave(1),:);
        w2t = Weights(ClassEstTemp == dirsave(2),:);
        scatter(d1t(:,otherdim),d1t(:,dimsave), w1t*scale,'ro')
        hold on
        scatter(d2t(:,otherdim),d2t(:,dimsave), w2t*scale,'bo')
        hold on
        SP=Tree.Layer{idt}.Threshold; %your point goes here 
        if(Tree.Layer{idt}.Dimension == 1)
            line([SP SP],get(hax,'YLim'),'Color',[0 0 1])
        else
            line(get(hax,'XLim'),[SP SP],'Color',[0 0 1])
        end
        hold off
        str = sprintf('adaboost training round %d',idt);
        title(str);
        xlabel(otherdim);
        ylabel(dimsave);
        pause(0.5)
        
    end
end

end
