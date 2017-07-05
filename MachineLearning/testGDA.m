%Test GDA


%create some fake data
u1 = [2,3];
u2 = [2,5];
c1 = [0.5,0.5];
c2 = [0.5,0.5];
nn = 100;

X0 = [c1(1)*randn(nn,1) + u1(1), c1(2)*randn(nn,1) + u1(2)];
X1 = [c2(1)*randn(nn,1) + u2(1), c2(2)*randn(nn,1) + u2(2)];

figure
scatter(X0(:,1),X0(:,2),'r.');
hold on
scatter(X1(:,1),X1(:,2),'b.');
hold off

%Get working data
X = [X0;X1];
Y = [ones(nn,1)*0; ones(nn,1)*1];




%Run GDA
[m,n] = size(X);
Y_row = repmat(Y,1,2);

phi = (1/m)*sum(Y==1);

mu_0 = sum( X.*(Y_row == 0) )/sum(Y==0); %find the rough mean of class 0
mu_1 = sum( X.*(Y_row == 1) )/sum(Y==1); %find the rough mean of class 1


%get relevant mean for each row
mu_y_0 = repmat(mu_0,m,1);
mu_y_1 = repmat(mu_1,m,1);
mu = mu_y_0.*(Y_row == 0) + mu_y_1.*(Y_row == 1);

sigma = (1/m).*( ( X - mu)'*( X - mu) );
sigin = inv(sigma);
thresh = 1/max(max(sigma));

%classify
stash = [];
for ii = 1:m
    P_0c = (1/ ( (2*pi)^(n/2) *sqrt(norm(sigma)) ) )* exp(-(1/2)*(X(ii,:) - mu_0)*sigin*(X(ii,:) - mu_0)' );
    P_1c = (1/ ( (2*pi)^(n/2) *sqrt(norm(sigma)) ) )* exp(-(1/2)*(X(ii,:) - mu_1)*sigin*(X(ii,:) - mu_1)' );
    
    if(max(P_0c/P_1c,P_1c/P_0c) > thresh)
        if(P_0c > P_1c)
            stash = [stash; Y(ii), 0, P_0c, P_1c, P_0c/P_1c];
        else
            stash = [stash; Y(ii), 1, P_0c, P_1c, P_1c/P_0c];
        end
    else
        %too similar, class unkown
        stash = [stash; Y(ii), -1, P_0c, P_1c, max(P_0c/P_1c,P_1c/P_0c)];
    end
end

subind = find(stash(:,2) ~= -1);
substash = stash(subind,:);
class = stash(:,1) == stash(:,2);
trueclass = substash(:,1) == substash(:,2);
disp('accuracy:')
acc = sum(class)/length(class)
disp('true accuracy:')
tacc = sum(trueclass)/length(trueclass)
disp('unkowns percentage:')
uacc = 1- (length(substash)/length(stash))

%For plotting
sigin = inv(sigma);
Xlims = [min(X),max(X)];
[X1p,X2p] = meshgrid(linspace(Xlims(1),Xlims(3),25)',linspace(Xlims(2),Xlims(4),25)');
Xp = [X1p(:) X2p(:)];

p = [];
for ii = 1:length(Xp)
    P_0temp = (1/ ( (2*pi)^(n/2) *sqrt(norm(sigma)) ) )* exp(-(1/2)*(Xp(ii,:) - mu_0)*sigin*(Xp(ii,:) - mu_0)' );
    P_1temp = (1/ ( (2*pi)^(n/2) *sqrt(norm(sigma)) ) )* exp(-(1/2)*(Xp(ii,:) - mu_1)*sigin*(Xp(ii,:) - mu_1)' );
    
    p = [p; P_0temp, P_1temp];
end

f1 = scatteredInterpolant(Xp(:,1),Xp(:,2),p(:,1));
P_0 = f1(X1p,X2p);
f1 = scatteredInterpolant(Xp(:,1),Xp(:,2),p(:,2));
P_1 = f1(X1p,X2p);
    
    
figure
scatter(X0(:,1),X0(:,2),'b.');
hold on
scatter(X1(:,1),X1(:,2),'r.');
hold on
handle0 = mesh(X1p,X2p,P_0);
hold on
handle1 = mesh(X1p,X2p,P_1);

hold off
set(handle0,'facecolor','none')
set(handle1,'facecolor','none')
colormap(cool);
colorbar;
    