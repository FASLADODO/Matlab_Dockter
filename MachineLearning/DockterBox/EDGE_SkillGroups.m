%% Group iterations by skill level 
%  Add k-fold partitions for k-fold crossvalidation later on.  

%DataKey.grp = grp;          % struct, use like this: grp.all{tsk}.Idx{g.Exp}
                            % or grp.prtn{tsk}.trn{g.Exp}{p}, & .eval{g.Exp}{p}
                            % .prtn sets up .prtnK partitions for
                            % k-fold cross validation.  For a skill
                            % group, 
                            % .eval{g}{p} is the left out set from the
                            % remaining .trn{g}{p}; e.g. .eval{g} is
                            % {[a], [b],[c]}, .trn{g} is
                            % {[b,c],[a,c],[a,b]}
%DataKey.grpFLS = grpFLS;    % same as .grp, but chosen based on FLS scores, 
                            % not demographics and experience
                            % Group names:
                            % g.Exp = 1, etc.  see above. 
if( ~exist('DataKey'))
    fprintf('Loading EdgeDataKey.mat ...');
    load EdgeDataKey.mat    
    fprintf(' DONE.');
end
load ../../raw-data/EDGE/scripts/vidReview.mat      
mrk={'rx', 'g.', 'bo'};                            
PLOT_GRP=1;

%% set & plot pOSATS and FLS thresholds for groundtruth expert group:            
pOSATS.thQ=3;    % Motion Qualtity score Q must be >= this AND  
pOSATS.thBi=3; % Bimanualtiy score must be >= this 
                 % for entrance to grndTruthExperts group (gtExp)

colBi=cellSearch('mnAll-Bi',vidRevLabels);
colQ=cellSearch('mnAll-Q',vidRevLabels);
colFLS=cellSearch('FLS',vidRevLabels);

if PLOT_GRP; figure(654);clf; end
for tsk=1:length(DataKey.Tasks)
    if PLOT_GRP;
    figure(654);
    plot(vidRev{tsk}.scores(:,colQ),vidRev{tsk}.scores(:,colBi),mrk{tsk});
    hold on;    
    end
    
    % figure(655)
    % subplot(1,3,tsk)
    % plot(vidRev{tsk}.scores(:,colQ),vidRev{tsk}.scores(:,colBi),mrk{tsk});
    % xlabel('Motion Quality');ylabel('Bimanuality');
    % title(DataKey.Tasks{tsk});
    % axis( [1 5 1 5])
    % axis square
    
    % find and save logs above p-OSATS thresholds
    gtExp{tsk}.logIdx=[];
    gtExp{tsk}.vidRevScoresIdx=[];
    gtExp{tsk}.vidRevScores=[];
    gtExp{tsk}.subjID=[];
    gtExp{tsk}.siteID=[];
    gtExp{tsk}.EXCLUDED.logIdx=[];
    gtExp{tsk}.EXCLUDED.subjID=[];
    gtExp{tsk}.EXCLUDED.siteID={};
    gtExp{tsk}.EXCLUDED.vidRevScoresIdx=[];
    gtExp{tsk}.EXCLUDED.vidRevScores=[];
    for i=1:size(vidRev{tsk}.scores,1);
        
        % if both scores above thresholds, save info
        if( vidRev{tsk}.scores(i,colQ)  >= pOSATS.thQ && ...
            vidRev{tsk}.scores(i,colBi) >= pOSATS.thBi ) 
        
            gtExp{tsk}.logIdx(end+1) = vidRev{tsk}.Idx(i);
            gtExp{tsk}.subjID(end+1) = DataKey.content{...
                vidRev{tsk}.Idx(i), DataKey.lookupCol('SubjID')};
            gtExp{tsk}.siteID{end+1} = DataKey.content{...
                vidRev{tsk}.Idx(i), DataKey.lookupCol('siteID')};
            gtExp{tsk}.vidRevScoresIdx(end+1)=[i];
            gtExp{tsk}.vidRevScores(end+1) =  vidRev{tsk}.scores(...
                i,  colFLS   );

        % otherwise, put it in .EXCLUDED
        else
            gtExp{tsk}.EXCLUDED.logIdx(end+1) = vidRev{tsk}.Idx(i);
            gtExp{tsk}.EXCLUDED.subjID(end+1) = DataKey.content{...
                vidRev{tsk}.Idx(i), DataKey.lookupCol('SubjID')};
            gtExp{tsk}.EXCLUDED.siteID{end+1} = DataKey.content{...
                vidRev{tsk}.Idx(i), DataKey.lookupCol('siteID')};
            gtExp{tsk}.EXCLUDED.vidRevScoresIdx(end+1)=[i];
            gtExp{tsk}.EXCLUDED.vidRevScores(end+1) =  vidRev{tsk}.scores(...
                i,  colFLS   );
        end
    end

    % Compute GrndTruthExperts_FLSthreshold
    gtFLSthreshold(tsk) = min( gtExp{tsk}.vidRevScores );
end

if PLOT_GRP
    figure(654)
    xlabel('Motion Quality');ylabel('Bimanuality');
    title('p-OSATS Scores');
    axis( [1 5.1 1 5.1])
    rectangle('position',[pOSATS.thQ pOSATS.thBi 5-pOSATS.thQ 5-pOSATS.thBi],...
        'EdgeColor',[1 0 0],'LineStyle','--' )
    txt = sprintf(['Ground Truth Expert Group (- - -)\n' ...
        'Min FLS scores = [ %.2f, %.2f, %.2f ]'],...
        gtFLSthreshold );
    text(2.1,1.5,txt,'color', [.8 .1 .1],'fontsize',7);
    disp(txt);
    
    axis square
    legend(DataKey.Tasks,'location', 'northwest','fontsize',8)
    
    toTexFig(['pOSATSGrpSlection2'] , ...
        ['Ground Truth expert group determination based on p-OSATS scores.' ], 4,1)
end


%% Build up gtExp plus group (holds same logs as gtExp but includes any
%  other logs above gtFLSthreshold they did);
taskCol = intersect(DataKey.lookupCol('DataLog'), DataKey.lookupCol('taskID'));
flsCol = DataKey.lookupCol('FLS-score');
for tsk=1:length(DataKey.Tasks)

    % for each subject (there should be only 1 log per 1 subject), save it,
    % then find any other logs by that subject whose scores exceed
    % gtFLSthreshold and save them too. 
    gtExpPlus{tsk}.logIdx = [];   
    allTaskIdx{tsk} = cellSearch(DataKey.Tasks{tsk},...
            DataKey.content(:,taskCol)');
    
    for i=1:length(gtExp{tsk}.logIdx)
        
        % save the original (pOSATS scored) log
        gtExpPlus{tsk}.logIdx(end+1) =  gtExp{tsk}.logIdx(i);
        
        % find all other logs by this subject (for this task)
        allSiteLogs = cellSearch(cell2mat(gtExp{tsk}.siteID(i)),...
            DataKey.content(:,DataKey.lookupCol('siteID'))' );
        
        allSubjIDLogs = cellSearch(gtExp{tsk}.subjID(i),...
            DataKey.content(:,DataKey.lookupCol('subjID'))' );
        
        allSubjLogs = intersect(allTaskIdx{tsk}, ...
            intersect(allSiteLogs, allSubjIDLogs));
        
        % find all logs over gtFLSthreshold
        for jj=allSubjLogs
            if (DataKey.content{jj,flsCol} >= gtFLSthreshold)
                
                % save them (append)
                gtExpPlus{tsk}.logIdx = ...
                    [gtExpPlus{tsk}.logIdx, jj];

            end
        end
    end
end




%% Build the group index (g) and data structure (grp) variables. 
g.gtExp=1; % ground truth (gt) experts: Faculty surgeons + Fellows--who got 
        % above pOSATS thesholds and best FLS-scoring log 
g.gtExpPlus=2; % same as grndTruthExperts with any of those 
        % subjects' logs that scored above gtFLSthreshold<task> added
        % gtFLSthreshold is the min FLS score of gtExperts<task>
g.flsExp=3; % all logs above gtFLSthreshold<task>
%g.flsPrf=4;
g.flsInt=4; % middle set of FLS scores (prcntile)
g.flsNov=5; % lowest percntile of FLS scores
g.flsMaster=6; % Very highest percntile of FLS scores (probably won't use)
g.Unknown=7; % all other logs that don't fit into other categories.  
g.n     =g.Unknown; % total count of groups
% NOTE: the combined logs from flsExp, flsInt, flsNov and Unknown should
% have no overlap and their union should hold ALL EverythingOK logs for
% that task.

grpTh.gtExp = gtFLSthreshold;      % by definition 
grpTh.gtExpPlus = gtFLSthreshold;  % by definition
grpTh.flsNovPrctile =  15; %only lowest 15 percentile of scores are Nov group
grpTh.flsIntPrctile =  []; % should be totally determined by flsNovPrctile

% grp.id{g.gtExp} should give 'Ground Truth Experts'
grp.id ={'GT-Experts', ...
    'GT-ExpertsPlus',...
    'FLS-Experts',...   %    'FLS-Proficient',...
    'FLS-Intermediate',...
    'FLS-Novice', ...
    'FLS-Master', ...
    'Unknown'}; %ISIS techs are 'unknown' ?
grp.tag={'gtExp', ...
    'gtExp+', ...
    'flsExp',...        %    'flsPrf',...
    'flsInt',...
    'flsNov',...
    'flsMst', ...
    'Unk' };

grp.all={};  % for given task, all log Idx numbers of a skill group
             % index like this: grp.all{tsk}.Idx{g.Exp}
grp.prtn={}; % for given task; partitions of .all for k-fold cross validation
             % index like this: grp.prtn{tsk}.trnIdx{g.Exp}{k} and .evalIdx{g.Exp}{k}
             % e.g., in a b c prtn, .trnIdx wil have ( b, ac, bc), .evalIdx (c,b,a)
grp.prtnK = 3;  % the k in k-fold cross validation, using 3. 
rseed = 209280 ;% found randomly with: sum(100*clock)
rand('seed', rseed); % seed rand for consistent results, for random partitions 

for tsk=1:length(DataKey.Tasks)
    gtExpIdx=[];
    gtExpPlusIdx=[];
    flsExpIdx=[]; %PrfIdx=[];
    flsIntIdx=[];
    flsNovIdx=[];
    flsMstIdx=[];
    UnkIdx=[];
    
    % gtExp group
    gtExpIdx = gtExp{tsk}.logIdx;
    
    % gtExpPlus group
    gtExpPlusIdx = gtExpPlus{tsk}.logIdx;
    
    % g.flsNov 's min FLS threshold 
    flsScores = cell2mat(DataKey.content(DataKey.LogIdx{tsk},...
        DataKey.lookupCol('FLS-score')));
    grpTh.flsNovMax = prctile(flsScores,grpTh.flsNovPrctile);
    %Note: this should be about equal: grpTh.flsNovPrctile*.01 ==?
    % length(find(flsScores <=grpTh.flsNovMax))/length(flsScores)
    
    % g.flsInt 's thresholds
    % compute mid: 'middle score' midpoint between Nov max and gtExp Min,
    % that is, between flsNovThreshold and gtFLSexpert threshold)
    mid = (grpTh.gtExp(tsk)-grpTh.flsNovMax)./2 + grpTh.flsNovMax;
    % now compute 'range' about the middle score:
    % try the same size of group based on expert size (not so good, cutting
    % group's experts range is HUGE!
    %rng = max(flsScores(find(flsScores>grpTh.gtExp(tsk)))) - ...
    %    min(flsScores(find(flsScores>grpTh.gtExp(tsk)))) ;
    % INSTEAD, compute the midpoint and grab N subjects around it
    % flsScore closest to mid:
    midIt = find([flsScores-mid]==min(abs(flsScores-mid))); 
%     % number of subject in flsInt group
%     Nint=round(length(flsScores)/20); % THIS (20) was used in training       
%     %Nint=10; 
%     % find range of scores of Nint scorers about mid
%     % compute range 
%     [s si] = sort(flsScores); % THIS  was used in training   
%     rng=abs(diff(s(find(si==midIt)+[-Nint Nint]))); % THIS  was used in training   
%     grpTh.flsIntMax = mid + rng/2;% .075;           % THIS  was used in training   
%     grpTh.flsIntMin = mid - rng/2;%  .075;          % THIS  was used in training   
    
    % instead: use same number of Novices for Intermediates
    numNovices = length(find(flsScores <=grpTh.flsNovMax));
    Nint = numNovices;
    [s si] = sort(flsScores); % THIS (20) was used in training   
    rng=abs(diff(s(find(si==midIt)+[-ceil(Nint/2) floor(Nint/2)]))); 
    grpTh.flsIntMax = mid + rng/2;%   
    grpTh.flsIntMin = mid - rng/2;

    
    % all other FLS-score based groups...
    for i = DataKey.LogIdx{tsk}' %IdxToPlot{tsk}'
        % for entries with clean data (this means allOK is true as well)
        % save the index to a list based on skill level
        
        flsScore = DataKey.content{i,DataKey.lookupCol('FLS-score')};
        
        % flsNov
        if( flsScore <= grpTh.flsNovMax)
            
            flsNovIdx = [flsNovIdx i];
        
        % flsExp    
        elseif (flsScore >= grpTh.gtExpPlus(tsk))
            
            flsExpIdx = [flsExpIdx i];
            
        % flsInt    
        elseif (flsScore <= grpTh.flsIntMax && flsScore >= grpTh.flsIntMin ) 
            
            flsIntIdx = [flsIntIdx i];
            
        else
            UnkIdx=[UnkIdx i];
            %warning('Unknown type of FLS score found ');
        end
        
    end
    
    % record the lists here.  Note: this should be indexable by 
    % grp.all{tsk}.Idx{g.Exp}() so (row) order is important
    grp.all{tsk}.Idx={gtExpIdx; ...
        gtExpPlusIdx;...
        flsExpIdx; ...%PrfIdx=[];
        flsIntIdx;...
        flsNovIdx;...
        flsMstIdx;...
        UnkIdx; };

    % partition logs randomly into k subgroups (exclusive bins) 
    % (skip unknown group)
    for gg=1:size(grp.all{tsk}.Idx, 1)-1
        
        n = length(grp.all{tsk}.Idx{gg});
        p = randperm(n);

        % split brackets
        first = round(linspace(1,n+1,grp.prtnK+1));
        last = first(2:end)-1;
        first(end)=[];

        % split indexes
        idx = arrayfun(@(i) p(first(i):last(i)), 1:grp.prtnK, 'uni', false);
        % split data
        partIdx = cellfun(@(r) grp.all{tsk}.Idx{gg}(r), idx', 'uni', false);

        % Check
        %partIdx{:}
        
        % save the .trn, and .eval,  evalIdx==the one left out in 
        % leave-one-out in k-fold cross validation; .trnIdx is what is
        % 'left
        % in'
        grp.prtn{tsk}.evalIdx{gg,:} = partIdx'; % like [{c}, {b}, {a}]
        for p=1:grp.prtnK
            grp.prtn{tsk}.trnIdx{gg}{p} = ... like [{a,b}, {a,c}, {b,c}]
                setdiff(  grp.all{tsk}.Idx{gg}, partIdx{p}); % note: order is important in setdiff
            %setxor( partIdx{p}', grp.all{tsk}.Idx{gg})
            
            % if the partitions overlap or don't include all members,
            % complain...
            if(~isempty(intersect(grp.prtn{tsk}.evalIdx{gg}{p},...
                    grp.prtn{tsk}.trnIdx{gg}{p})) || ~isempty( setdiff(...
                union(grp.prtn{tsk}.evalIdx{gg}{p},grp.prtn{tsk}.trnIdx{gg}{p}),...
                grp.all{tsk}.Idx{gg})))
                error('Failed to setup k-fold partitioning correctly')
            end
        end
    end
end



%% Plot Groups here.
% nevermind,  Made plots down below (after FLS group plots)




%% Group iterations by FLS skill level and create  k-fold partitions 
%  for k-fold cross validation later on. 

%% FLS Group Boundaries (based on FLS scores, the ultra top, top, middle and 
% bottom scorers
flsGrpTh = [90 100; 80 100; 40 60; 0.01 20; ]; % percent of scorers: 0-20:Nov,  
                                % 40-60:Int,  
%g.Exp=1; % FLS Experts group-- very top, 90-100% of scorers
%g.Prf=2; % FLS Proficient Group. -- top, 80-100% of scorers
%g.Int=3; % FLS Intermediate Group -- mid, 40-60% of scorers
%g.Nov=4; % FLS Novices -- bottom, 0-20% of scorers
%g.FLSmst=5; % Masters -- These are the superstars at each site, HUGE amount of experience
grpFLS.id ={'FLS-Expert', 'FLS-Proficient', 'FLS-Intermediate', 'FLS-Novice'};%,'Unknown'}; %ISIS techs are 'unknown' classification
grpFLS.tag={'ExpFLS',     'PrfExp',         'IntFLS',           'NovFLS' };%, 'Unk' };
grpFLS.all={};  % for given task, compat set of all Idx numbers of a skill group
             % index like this: grpFLS.all{tsk}.Idx{g.Exp}
grpFLS.prtn={}; % for given task; partitions of .all for k-fold cross validation
             % index like this: grpFLS.prtn{tsk}.trnIdx{g.Exp}{k} and .evalIdx{g.Exp}{k}
             % e.g., in a b c prtn, .trnIdx wil have ( b, ac, bc), .evalIdx
             % (c,b,a)
grpFLS.prtnK=grp.prtnK;  % k in k-fold cross validation, use same as .grp
flsScoreCol = intersect(DataKey.lookupCol('Summary Metrics'), ...
    DataKey.lookupCol('FLS-Score'));

colordef white
% first extract FLS score for each task, each log
for tsk=1:length(DataKey.Tasks)   
    FLSscores{tsk}.all = cell2mat(...
        DataKey.content(DataKey.LogIdx{tsk} , flsScoreCol) );    
    
    N=length(FLSscores{tsk}.all);
    [s sortIdx{tsk}]=sort(FLSscores{tsk}.all);
    
    for gg=1:size(flsGrpTh,1)
        idxRng = round(prctile(1:N, flsGrpTh(gg,:)));
        scoresRng= s(idxRng(1):idxRng(end));
        idxGrpAll=[];
        idxAllScore=[];
        for i = DataKey.LogIdx{tsk}'    
            if(any(cell2mat(DataKey.content(i , flsScoreCol))==...
                   scoresRng))
                idxGrpAll(end+1) = i;                 
            end
            idxAllScore(end+1) = i;
        end
        grpFLS.all{tsk}.Idx{gg,:} = idxGrpAll;    % holds all Idx's (excel row#'s) for a group
        %grpFLS.all{tsk}.IdxScore{gg,:} = idxGrpAllScore; % holds corresponding score for those Idx's
    end
    grpFLS.scores{tsk}=FLSscores{tsk}.all'; % Dense array of scores, no indexing here.
    FLSscores{tsk}.all= FLSscores{tsk}.all';
    grpFLS.scoresIdx{tsk}= idxAllScore; % mapping to Idx of dense score array.  
    FLSscoresIdx{tsk} = idxAllScore;
    
    if PLOT_GRP
        fontsize=8;
        figure(tsk+50)
        subplot(3,1,1:2)
        nbins=50;
        [n x]=hist(grpFLS.scores{tsk},nbins);
        hist(grpFLS.scores{tsk},50);
        title(['\bfFLS Scores - ' DataKey.Tasks{tsk} ...
            ' (N=' num2str(N) ')' ], 'fontsize', fontsize)
        ylabel('Count', 'fontsize', fontsize)
        axis([-1 1.2, 0 18])
        %set(gca, 'Units','normalized', 'Position',[0.15 0.2 0.75 0.7]);
        op = get(gca, 'outerposition')
        op(1)=0; op(3)=1;
        set(gca,'OuterPosition' , op)
        
        
        subplot(3,1,3)
        % e=[-1000:.01:1];
        % n = histc(FLSscores{tsk}.all,e);
        %plot(e,cumsum(n)/max(sum(n)),'.');
        %plot(s,cumsum(s)/sum(s),'.')
        plot(x, cumsum(n)./sum(1),'.')
        ylabel('Cumulative Count', 'fontsize', fontsize)
        xlabel('FLS Score', 'fontsize', fontsize)
        axis([-1 1.2, 0 200])
        %set(gca, 'Units','normalized', 'Position',[0.15 0.2 0.75 0.7]);
        op = get(gca, 'outerposition')
        op(1)=0; op(3)=1;
        set(gca,'OuterPosition' , op)
        
        pos = get(gcf,'position');
        pos(3)=3/4*pos(4);
        set(gcf,'position',    pos)
        %toTexFig(['FlsScoreHistogram' DataKey.Tasks{tsk}] , ...
        %         ['FLS Scores,  ' DataKey.Tasks{tsk} ' task'], 2,1);
        w=1.8; h=2.5;
        
        set(gcf, 'units', 'inches')
        set(gcf, 'position', [ 0 0 w h])
        toTexFigPageSize(['FlsScoreHistogram' DataKey.Tasks{tsk}] , ...
             ['FLS Scores,  ' DataKey.Tasks{tsk} ' task'], [w h]);
    end
    
    % partition FLS group logs randomly into k subgroups (exclusive bins) 
    % (skip unknown group)
    for gg=1:size(flsGrpTh,1)
        
        n = length(grpFLS.all{tsk}.Idx{gg});
        p = randperm(n);

        % split brackets
        first = round(linspace(1,n+1,grpFLS.prtnK+1));
        last = first(2:end)-1;
        first(end)=[];

        % split indexes
        idx = arrayfun(@(i) p(first(i):last(i)), 1:grpFLS.prtnK, 'uni', false);
        % split data
        partIdx = cellfun(@(r) grpFLS.all{tsk}.Idx{gg}(r), idx', 'uni', false);

        % Check
        %partIdx{:}
        
        % save the .trn, and .eval,  evalIdx==the one left out in 
        % leave-one-out in k-fold cross validation; .trnIdx is what is
        % 'left
        % in'
        grpFLS.prtn{tsk}.evalIdx{gg,:} = partIdx'; % like [{c}, {b}, {a}]
        for p=1:grpFLS.prtnK
            grpFLS.prtn{tsk}.trnIdx{gg}{p} = ... like [{a,b}, {a,c}, {b,c}]
                setdiff(  grpFLS.all{tsk}.Idx{gg}, partIdx{p}); % note: order is important in setdiff
            %setxor( partIdx{p}', grpFLS.all{tsk}.Idx{gg})
            
            % if the partitions overlap or don't include all members,
            % complain...
            if(~isempty(intersect(grpFLS.prtn{tsk}.evalIdx{gg}{p},...
                    grpFLS.prtn{tsk}.trnIdx{gg}{p})) || ~isempty( setdiff(...
                union(grpFLS.prtn{tsk}.evalIdx{gg}{p},grpFLS.prtn{tsk}.trnIdx{gg}{p}),...
                grpFLS.all{tsk}.Idx{gg})))
                error('Failed to setup k-fold partitioning correctly')
            end
        end
    end
    
end
grpFLS.sortIdx=sortIdx;
FLSscoresSortIdx=sortIdx;

%% Plot sorted scores for all FLS tasks (FLS GROUPS)
mrk={'r.', 'go', 'mo', 'k.', 'co'};
mrkTsk = {'k--', 'k-', 'k:'};
if PLOT_GRP
    grpToPlot=fliplr([ g.flsMaster g.flsExp  g.flsInt g.flsNov]);
    figure
    for tsk=1:length(DataKey.Tasks)
        plot(  grpFLS.scores{tsk}(grpFLS.sortIdx{tsk}), mrkTsk{tsk});hold on;
    end
    legend(DataKey.Tasks,'location', 'best')
    
    for tsk=1:length(DataKey.Tasks)
        for gg=4:-1:1%length( grpFLS.all{tsk}.Idx)
            ii=find(ismember(grpFLS.scoresIdx{tsk}',...
                grpFLS.all{tsk}.Idx{gg,:}));
            iiSort=find(ismember(grpFLS.sortIdx{tsk},ii));
            plot( iiSort,...
                sort(grpFLS.scores{tsk}(ii)),mrk{gg});
        end
        if(tsk==1)
            h=legend({DataKey.Tasks{:},grpFLS.id{4:-1:1}},'location', 'southeast');
            set(h, 'fontsize', 7)
        end
        %    for gg=1:length( grpFLS.all{tsk}.Idx)
        %        plot( grpFLS.all{tsk}.Idx{gg},mrk{gg})
        %    end
    end
    xlabel('Iteration Sorted by FLS Score');
    ylabel('FLS Score');
    title(['FLS Percentile-Based Group Selection']);
    
    
    set(gcf,'position',    get(gcf,'position')*diag([ 1 1 1 1]))
    
    toTexFig(['FlsGrpSlection2' ] , ...
        ['FLS Score-based group selection' ], 4,1)
    
end
 %%

 
 
 %% Plot sorted scores for all tasks (Ground Truth Group)
 mrk={'r.', 'bo', 'go', 'mo', 'k.'};
 mrkTsk = {'k--', 'k-', 'k:'};
 %grpToPlot=[g.flsNov g.flsInt g.flsExp  g.gtExp];
 grpToPlot=fliplr([ g.gtExp g.gtExpPlus g.flsExp  g.flsInt g.flsNov]);
 if PLOT_GRP
     figure
     for tsk=1:length(DataKey.Tasks)
         plot(  FLSscores{tsk}.all(FLSscoresSortIdx{tsk}), mrkTsk{tsk});hold on;
     end
     legend(DataKey.Tasks,'location', 'best')
     for tsk=1:length(DataKey.Tasks)
         for gg=grpToPlot%4:-1:1%length( grpFLS.all{tsk}.Idx)
             ii=find(ismember(FLSscoresIdx{tsk}',...
                 grp.all{tsk}.Idx{gg,:}));
             %iiSort=find(ismember(grp.sortIdx{tsk},ii));
             iiSort=find(ismember(FLSscoresSortIdx{tsk},ii));
             plot( iiSort,...
                 sort(FLSscores{tsk}.all(ii)),mrk{gg})
         end
         if(tsk==1)
             h=legend({DataKey.Tasks{:},grp.id{grpToPlot}},'location', 'southeast');
             set(h, 'fontsize', 7)
         end
         %    for gg=1:length( grpFLS.all{tsk}.Idx)
         %        plot( grpFLS.all{tsk}.Idx{gg},mrk{gg})
         %    end
     end
     xlabel('Iteration Sorted by FLS Score');
     ylabel('FLS Score');
     title(['Ground Truth-Based Data Selection (per iteration)']);
     
     
     %set(gcf,'position',    get(gcf,'position')*diag([ 1 1 1 1]))
     
     toTexFig(['GroundTruthGrpSlection2'] , ...
         ['FLS Score-based group selection' ], 4,1)
     drawnow;
 end


%% Make table of Ground Truth Expert and Group selection

% fprintf('\n\n')
% fprintf('\n\n\\toprule\n')
% for tsk=1:length(DataGlb.Tasks)
%    fprintf('  & \\bfseries %s  ',DataGlb.Tasks{tsk}); 
% end
% fprintf(' \\\\ \n');
% fprintf('\\midrule\n');
% 
% for i=1:length(toPlot)
%     fprintf('\\multirow{2}{*}{%6s} & $r$ \t& %.2f (%.2f) & %.2f (%.2f) & %.2f (%.2f)\\\\ \n',...
%         xLabels{i},pearsonR(i,:))
%     fprintf('                            & $\\rho$ \t& %.2f (%.2f) & %.2f (%.2f) & %.2f (%.2f) \\\\ \n',...
%         spearmanR(i,:))
%     fprintf('\\cmidrule(r){2-5}\n');
% end
% fprintf('\\bottomrule\n');
% fprintf('\n\n');
% 
% 
% 
% 
% for i=1:length(toPlot)
%     fprintf('%s ',DataGlb.content{2,toPlot(i)});
%     for tsk=1:length(DataGlb.Tasks)
%         if tsk==0
%             fprintf('$%0.2f(p=%.2f)$ ',spearmanR(i,tsk),spearmanRp(i,tsk));
%         else
%             fprintf(' & $%0.2f(%.2f)$',spearmanR(i,tsk),spearmanRp(i,tsk));
%         end        
%     end
%     fprintf(' \\\\\n')
% end
% fprintf('\\bottomrule\n')