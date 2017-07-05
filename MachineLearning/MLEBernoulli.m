%% MLE bernoulli

%Simulate some bernouli data
N = ones(1000,1);
P = 0.3 %acutal probability

R = binornd(N,P);

%Actuall probabilities
P0 = sum(R == 0) / length(N)
P1 = sum(R == 1) / length(N)


%Likelihood
%http://www.colorado.edu/economics/morey/7818/estimation/maxlik/maxlik.pdf

psam = 0:0.01:1;
NS = length(N);
xbar = mean(R);

LnL = NS*xbar.*log(psam) + NS*(1-xbar).*log(1-psam);

[Lml, idx] = max(LnL);

pml = psam(idx)

figure
plot(psam,LnL)
hold on
scatter(pml,Lml,'g*')
hold off