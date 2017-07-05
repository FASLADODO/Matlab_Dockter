% Entropy vs variance

nC = 6; %number of classes
Probs = [];

% Even
Probs(end + 1,:) = ones(1,nC)/nC;

%zeros + ones in places
for jj = 1:nC
    Probs(end + 1,:) = zeros(1,nC);
    Probs(end,jj) = 1;
end

%1/n and (n-1)/n
for ii = 1:nC
    for jj = ii:nC
        if( ii ~= jj)
            Probs(end + 1,:) = zeros(1,nC);
            Probs(end,ii) = (nC-1)/nC;
            Probs(end,jj) = 1/nC;
        end
    end
end

%2/n and (n-2)/n
for ii = 1:nC
    for jj = ii:nC
        if( ii ~= jj)
            Probs(end + 1,:) = zeros(1,nC);
            Probs(end,ii) = (nC-2)/nC;
            Probs(end,jj) = 2/nC;
        end
    end
end


%halves
for ii = 1:nC
    for jj = ii:nC
        if( ii ~= jj)
            Probs(end + 1,:) = zeros(1,nC);
            Probs(end,ii) = 1/2;
            Probs(end,jj) = 1/2;
        end
    end
end

%thirds
for ii = 1:nC
    for jj = ii:nC
        if( ii ~= jj)
            for kk = jj:nC
                if(jj ~= kk)
                    Probs(end + 1,:) = zeros(1,nC);
                    Probs(end,ii) = 1/3;
                    Probs(end,jj) = 1/3;
                    Probs(end,kk) = 1/3;
                end
            end
        end
    end
end

%fourths
for ii = 1:nC
    for jj = ii:nC
        if( ii ~= jj)
            for kk = jj:nC
                if(jj ~= kk)
                    for ll = kk:nC
                        if(kk ~= ll)
                            Probs(end + 1,:) = zeros(1,nC);
                            Probs(end,ii) = 1/4;
                            Probs(end,jj) = 1/4;
                            Probs(end,kk) = 1/4;
                            Probs(end,ll) = 1/4;
                        end
                    end
                end
            end
        end
    end
end

% Alocate metric vectors
CombineVar = [];
CombineEnt = [];
EntropyVec = [];
VarVec = [];
[Row,Col] = size(Probs);


% Get entropy
for ii = 1:Row
   EntropyVec(ii) = EntropyCalc(Probs(ii,:)); 
end

% Get Variance
for ii = 1:Row
   VarVec(ii) = Variance(Probs(ii,:)); 
end

%Combined probabilities 
% for ii = 1:Row
%     for jj = ii:Row
%         if(ii ~= jj)
%             CombineEnt(end+1) = IntraClassEntropy(Probs(ii,:) , Probs(jj,:));
%         end
%     end
% end

%% Plot sample diff

% sample = [1,53,81];%generic
sample = [53,61,72];%movement

count = 1;
for ii = 1:length(sample)
    for jj = ii:length(sample)
        kk = sample(ii);
        ll = sample(jj);
        CombineEnt(count) = IntraClassEntropy(Probs(kk,:) , Probs(ll,:));
        CombineVar(count) = IntraClassVar(Probs(kk,:) , Probs(ll,:));
        
        figure(count)
        subplot(1,4,1);
        bar(Probs(kk,:),0.9,'r')
        hold on 
        bar(Probs(ll,:),0.5,'b')
        hold off
        legend('-A','-B')
        title('PDFs')
        xlabel('Outcome')
        ylabel('Probabilities')

        subplot(1,4,2);
        bar([EntropyVec(kk),0],0.9,'r')
        hold on 
        bar([0,EntropyVec(ll)],0.5,'b')
        hold off
        set(gca,'XTickLabel',{'Class A', 'Class B'})
        legend('-A','-B')
        title('Entropy')
        ylabel('Entropy')
        
        subplot(1,4,3);
        bar([VarVec(kk),0],0.9,'r')
        hold on 
        bar([0,VarVec(ll)],0.5,'b')
        hold off
        set(gca,'XTickLabel',{'Class A', 'Class B'})
        legend('-A','-B')
        title('Variance')
        ylabel('Variance')
        
        subplot(1,4,4);
        bar([CombineEnt(count),0],'g')
        hold on
        bar([0,CombineVar(count)],'c')
        hold off
        set(gca,'XTickLabel',{'Entropy', 'Variance'})
        title('Combined')
        ylabel('Var/Ent')

        str = sprintf('Probability %i vs %i',kk,ll);
        suptitle(str)
        count = count + 1;
    end
end

%% Plot thresh

count = 1;
for kk = 1:Row
    minner = 10000;
    mindex = 0;
    maxxer = 0;
    maxdex = 0;
    for ll = kk:Row
        IntraEnt = IntraClassEntropy(Probs(kk,:) , Probs(ll,:));
        IntraVar = IntraClassVar(Probs(kk,:) , Probs(ll,:));
        
        checkRatio = 0;
        if(EntropyVec(kk) ~= 0 && EntropyVec(ll) ~= 0 && IntraEnt ~= 0)
            checkRatio = (IntraEnt/(EntropyVec(kk)*EntropyVec(ll)));
            if(checkRatio > maxxer)
               maxxer = checkRatio; 
               maxdex = ll;
            end
            if(checkRatio < minner)
               minner = checkRatio; 
               mindex = ll;
            end
        end
        
        if( checkRatio ~=0 )
%             CombineEnt(count) = IntraEnt;
%             CombineVar(count) = IntraVar;
% 
%             figure(count)
%             subplot(1,4,1);
%             bar(Probs(kk,:),0.9,'r')
%             hold on 
%             bar(Probs(ll,:),0.5,'b')
%             hold off
%             legend('-A','-B')
%             title('PDFs')
%             xlabel('Outcome')
%             ylabel('Probabilities')
% 
%             subplot(1,4,2);
%             bar([EntropyVec(kk),0],0.9,'r')
%             hold on 
%             bar([0,EntropyVec(ll)],0.5,'b')
%             hold off
%             set(gca,'XTickLabel',{'Class A', 'Class B'})
%             legend('-A','-B')
%             title('Entropy')
%             ylabel('Entropy')
% 
%             subplot(1,4,3);
%             bar([VarVec(kk),0],0.9,'r')
%             hold on 
%             bar([0,VarVec(ll)],0.5,'b')
%             hold off
%             set(gca,'XTickLabel',{'Class A', 'Class B'})
%             legend('-A','-B')
%             title('Variance')
%             ylabel('Variance')
% 
%             subplot(1,4,4);
%             bar([CombineEnt(count),0],'g')
%             hold on
%             bar([0,CombineVar(count)],'c')
%             hold off
%             set(gca,'XTickLabel',{'Entropy', 'Variance'})
%             title('Combined')
%             ylabel('Var/Ent')
% 
%             str = sprintf('Probability %i vs %i',kk,ll);
%             suptitle(str)
%             count = count + 1;
        end
    end
    resultsArray(kk,:) = [maxxer,maxdex,minner,mindex];
end



