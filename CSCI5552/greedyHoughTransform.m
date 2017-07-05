clear all

%% greedy hough transform

nn = 10;

x1 = [linspace(2,2.5,nn)'];
x2 = [linspace(1,10,nn)'] + randn(nn,1)*0.2;

y1 = [linspace(1,10,nn)'] + randn(nn,1)*0.2;
y2 = [linspace(1,1.2,nn)'];



figure
scatter(x1(:,1),y1,'r*')
hold on 
scatter(x2(:,1),y2,'b+')
hold off


%% hough

%data = x, 1, y

data = [x1,y1; x2, y2];
numsensors = size(data,1);

threshvote = 1;

thetarng = linspace(deg2rad(-90),deg2rad(90),10);
mindist = -max(abs(NormRowWise([data(:,1),data(:,2)])));
maxdist = max(abs(NormRowWise([data(:,1),data(:,2)])));
distrng = linspace(mindist,maxdist,10);


accumulator = zeros(length(thetarng),length(distrng));
bookeeper = zeros(length(thetarng),length(distrng),numsensors);

%Go through each point and add votes to accumulator array
tic
ids = 1;
store = [];
for ii = 1:numsensors
   for jj = 1:length(thetarng) %test each angle
       rt = data(ii,1)*cos(thetarng(jj)) + data(ii,2)*sin(thetarng(jj)); %compute r
       store(ids) = rt; %for checking later
       ids = ids + 1;
       if(rt > maxdist)
          rt = maxdist - 1;
          disp('big')
       elseif(rt < mindist)
          rt = mindist + 1; 
          disp('little')
       end
       %find which distance bin to put in
       id = find(distrng > rt,1,'first');
       bookeeper(jj,id,ii) = 1;
   end
end
toc

%plot initial vote surface
accumulatorcpy = sum(bookeeper,3);
figure(10);
handle = surf(accumulatorcpy);
xlabel('dist')
ylabel('theta')
zlabel('votes')
ax = gca;

ax.XTick = 1:length(distrng);
ax.YTick = 1:length(thetarng);
ax.XTickLabel = round(distrng);
ax.YTickLabel = thetarng;

title('accumulator array 1')

%%

lines = [];
ids = 1;

tic
while(1)
    if(1)
        %plot initial vote surface
        accumulator = sum(bookeeper,3);
        figure
        handle = surf(accumulator);
        xlabel('dist')
        ylabel('theta')
        zlabel('votes')
        ax = gca;
        ax.XTickLabel = distrng;
        ax.YTickLabel = thetarng;
        title('accumulator array')
    end
    
    %find the max from bookeeper accumulator
    accumulator = sum(bookeeper,3);
    topvote = max(max(accumulator))
    if(topvote <= threshvote)
        break;
    end
    
    %find the element in accumulator with the highest votes
    mid = find(accumulator == topvote, 1,'first');
    [rm,cm] = ind2sub(size(accumulator),mid);
    
    %get the vote stack from 3D matrix
    arr = bookeeper(rm,cm,:);
    idarr = permute(arr,[3,2,1]) %conver to 1D array
    
    %get data for points that voted for this line
    ppd = [data(idarr == 1,1), ones(length(data(idarr == 1,1)), 1)];
    ppy = data(idarr == 1,2);
    
    %compute actual line from these points
    params = pinv(ppd)*ppy; 
    
    %store it
    lines(:,ids) = params; %m, b
    ids = ids + 1;
    
    % now remove those indices from the bookeeper
    checker = repmat(arr,length(thetarng),length(distrng),1);
    bookeeper(bookeeper == checker) = 0;
    
    
end
toc


%% plotting


xp = [linspace(1,8,nn)', ones(nn,1)];

figure
for ii = 1:size(lines,2)
    yp = xp*lines(:,ii)
    plot(xp(:,1),yp)
    hold on
end

scatter(x1(:,1),y1,'r*')
hold on 
scatter(x2(:,1),y2,'b+')
hold off


