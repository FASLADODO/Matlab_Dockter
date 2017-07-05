%testCCA

%create some fake data
mu1 = [1,2];
sigma1 = [2,1; 1,2];

mu2 = [2,3];
sigma2 = [3,1; 1,3];

nn = 100;

X0 = mvnrnd(mu1,sigma1,nn);
X1 = mvnrnd(mu2,sigma2,nn);

[A,B,r,U,V]= canoncorr(X0,X1)

figure
plot(U(:,1),V(:,2),'.')


%Get working data
X = [X0;X1];
Labels = [ones(nn,1)*-1; ones(nn,1)*1]; %(-1/1)


figure
gscatter(X(:,1),X(:,2),Labels);


%% feature selection

mdl = fscnca(X,y,'Solver','sgd','Verbose',1);

%%

[NN,SS] = size(X);
if(SS > NN)
    X = X'; 
    [NN,SS] = size(X);
end


cslist = unique(Labels);

%get data
for ii = 1:length(cslist)
    idx = Labels == cslist(ii);
    X_Class{ii} = X(idx,:);
end

%get cross correlations
for ii = 1:length(cslist)
    Varz{ii} = cov(X_Class{ii});
    for jj = 1:length(cslist)
        Covz{ii,jj} = cov(X_Class{ii},X_Class{jj});
    end
end

Covz{1,2} ./ sqrt(Varz{1}.*Varz{2})

