%% BISHOP Book


N = 10;
noiz = 0.2;

x = linspace(0,1,N)';

t = sin(2*pi.*x) + randn(N,1)*noiz;

scatter(x,t,'r.')

%% Try polynomial fitting for various orders and find the minimum error

M = 0:9;

%For finer plotting
NF = 100;
xf = linspace(0,1,NF)';

ttest = sin(2*pi.*xf) + randn(NF,1)*0.3;

for mm = 1:length(M)

    D = ones(N,1);
    Donline = ones(NF,1);
    for ii = 1:M(mm)
        D = [D, x.^ii];
        Donline = [Donline, xf.^ii];
    end
    params = pinv(D)*t;

    %Get out new t
    ft = Donline*params;

    sct = D*params;

    RMSError(mm) = ErrorRMS(ft,ttest);

end

figure
plot(RMSError)

%plot em
figure
scatter(x,t,'r.')
hold on
plot(xf,ft,'b-')
hold off

%% Probabilities

Pbox = [4/10; 6/10]; %[red; blue]

PJ = [1/4, 3/4;3/4,1/4] % [r; b] [a, o]

P_fruit = PJ*Pbox

%Bayes
P_cond = PJ.*[Pbox,Pbox] ./ [P_fruit,P_fruit]'

%% Bootstrap error bars

M = 4; %the minimum from the orders test for polynomial

BS = 20;
PP = 0.4;

%For finer plotting
NF = 100;
xf = linspace(0,1,NF)';

ttest = sin(2*pi.*xf) + randn(NF,1)*0.1;


for ii = 1:BS
    rng('shuffle')
    idxr = sort( randsample(NF,round(PP*NF)) );
    
    xon = xf(idxr,:);
    ton = ttest(idxr,:);
    
    Don = ones(length(xon),1);
    for kk = 1:M
        Don = [Don, xon.^kk];
    end
    params = pinv(Don)*ton;
    
    saveparams(:,ii) = params;
    
    figure
    plot(xon,ton)
end

figure
for jj = 1:M+1
scatter(ones(1,BS)*jj,saveparams(jj,:),'r.')
hold on
end
hold off

uncert = var(saveparams,[],2)
W = mean(saveparams,2)

figure
errorbar([1:(M+1)],W,uncert)

%% Entropy and stuff

P = ones(8,1)/8
H = Entropy(P)


P2 = [1/2,1/4,1/8,1/16,1/64,1/64,1/64,1/64]
H2 = Entropy(P2)

%% Prob 1.3


Pbox = [2/10; 2/10; 6/10]; %[red; blue green]

PJ = [3/10, 4/10, 3/10; 1/2, 1/2, 0/2; 3/10, 3/10, 4/10] % [r; b; g] [apple, orange, lime]

P_fruit = PJ'*Pbox

%Bayes
P_cond = PJ.*[Pbox,Pbox,Pbox] ./ [P_fruit,P_fruit,P_fruit]
