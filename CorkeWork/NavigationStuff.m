%% Navigation

[x,y] = meshgrid(0:0.1:100,0:0.1:100); 

s = sensorfield(x,y);

sim('sl_braitenberg')

figure(1)
h=surf(x,y,s);
set(h,'LineStyle','none')
view(2)

hold on

plot(xout(:,1),xout(:,2))
hold off


%% bugs

load map1

bug = Bug2(map);

bug.goal = [50,35]; %end point

bug.path([30,90]); %starting point

%% make maps

map = makemap(100);

bug = Bug2(map);

bug.goal = [30,90]; %end point

bug.path([90,10]); %starting point

%% D* navigation

load map1
goal = [50;30];
start = [20;10];

ds = Dstar(map);

ds.plan(goal);

ds.path(start);


%% Voronoi nav

load map1
goal = [50;30];
start = [20;10];

%find the free space
free = 1 - map;
free(1,:) = 0;
free(100,:) = 0;
free(:,1) = 0;
free(:,100) = 0;

skeleton = ithin(free);
imshow(skeleton)

%% probabalistic road map

load map1
goal = [50;30];
start = [20;10];

prm = PRM(map);

prm.plan();

prm.visualize();

prm.path(start,goal)

%% Rapidly-exploring Random Tree (RRT)
goal = [0,0,0];
start = [0,2,0];
veh = Vehicle([], 'stlim', 1.2);
rrt = RRT([], veh, 'goal', goal, 'range', 5);
rrt.plan()             % create navigation tree
rrt.path(start, goal)  % animate path from this start location
