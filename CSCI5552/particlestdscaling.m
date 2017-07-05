offset = 5;
alpha = 4;
beta = 3;

bestw = 0:0.01:1;

sigma = alpha.*exp(beta.*bestw) + offset;

figure
scatter(bestw,sigma)

% figure
% scatter(bestw,exp(bestw))