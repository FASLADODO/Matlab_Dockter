
% Grasp speed
myGroups = fliplr([g.gtExp g.flsInt g.flsNov ]);
%myGroups = g.gtExp ;
GraspData = [];

tskG = [1,2,3];

% get all grasp data from all grasps ever
for gg = 1:length(myGroups)
    for tg = tskG
        % Sum up all logs specific group (int, nov, exp)
        for ii = 1:length(DataGlb.grp.all{tg}.Idx{myGroups(gg)})

            i = DataGlb.grp.all{tg}.Idx{myGroups(gg)}(ii);

            %Get grasp data
            GR = DataGlb.dataLog{i}(:,G.QgR);
            GL = DataGlb.dataLog{i}(:,G.QgL);
            dGR = DataGlb.dataLog{i}(:,G.dQgR);
            dGL = DataGlb.dataLog{i}(:,G.dQgL);

            GraspData = [GraspData; dGR,dGL];
        end
    end
end

%get a histogram
X = abs(GraspData(:));
figure
histogram(X)


%Now sort to get 98 % of all grasp speeds
a = 0.98;
XS = sort(X);

nn = length(XS);

XP = XS(1:round(nn*a),:);

XP(end)