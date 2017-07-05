%k means from scratch

clear all


%create some fake data
u1 = [2,3];
u2 = [3,4];
u3 = [2,5];
u4 = [10,15];
c1 = [0.5,0.5];
c2 = [0.5,0.5];
c3 = [0.5,0.5];
c4 = [0.5,0.5];
nn = 100;

X0 = [c1(1)*randn(nn,1) + u1(1), c1(2)*randn(nn,1) + u1(2)];
X1 = [c2(1)*randn(nn,1) + u2(1), c2(2)*randn(nn,1) + u2(2)];
X2 = [c3(1)*randn(nn,1) + u3(1), c3(2)*randn(nn,1) + u3(2)];
X3 = [c4(1)*randn(nn,1) + u4(1), c4(2)*randn(nn,1) + u4(2)];


figure(1)
scatter(X0(:,1),X0(:,2),'r.');
hold on
scatter(X1(:,1),X1(:,2),'b.');
hold on
scatter(X2(:,1),X2(:,2),'c.');
hold on
scatter(X3(:,1),X3(:,2),'k.');
hold off

%Get working data
X = [X0;X1;X2;X3];
Y = [ones(nn,1)*1; ones(nn,1)*2; ones(nn,1)*3; ones(nn,1)*4];

%% weaponize it

km = 2;

[means,grpIndex] = kMeansIterative(X,km);

figure
for mm = 1:km
    groupdata = X(grpIndex == mm, :);

    scatter(groupdata(:,1),groupdata(:,2))
    hold on
    scatter(means(mm,1),means(mm,2),'g*')
    hold on
end
hold off

%%

[NN,SS] = size(X);
kmeans = 3;
maxiter = 1000;

%initialize means
rng('shuffle')
seeds = randi(NN, kmeans, 1);
means = X(seeds,:)

% figure
% scatter(X0(:,1),X0(:,2),'r.');
% hold on
% scatter(X1(:,1),X1(:,2),'b.');
% hold on
% scatter(means(:,1),means(:,2),'g*');
% hold off


clusters = [];

groups = ones(NN,1);
prevgroups = groups;

doneso = 0;
figure(2)
for tt = 1:maxiter

    %compute voronoi groups
    dist2means = zeros(NN,kmeans);

    for mm = 1:kmeans
        scale = X - repmat(means(mm,:),NN,1);
        dist = sqrt(sum(abs(scale).^2,2));
        
        dist2means(:,mm) = dist;
    end

    
    [~,groups] = min(dist2means,[],2);


    clf
    %get data groupings and compute new means
    for mm = 1:kmeans
        groupdata = X(groups == mm, :);
        
        scatter(groupdata(:,1),groupdata(:,2))
        hold on
        scatter(means(mm,1),means(mm,2),'g*')
        hold on
        
        means(mm,:) = mean(groupdata,1);
    end
    hold off
    w = waitforbuttonpress;
    if w == 0
        break;
    else
        
    end
    
    %check if we stopped changing
    if(isequal(groups,prevgroups) )
       doneso = 1;
       fprintf('convergence in %d steps',tt)
       break; 
    end
    prevgroups = groups;
end

means

if(doneso == 0)
   disp('failed to converge!') 
end
