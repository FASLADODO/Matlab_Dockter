%% Evaluate sample gausian signal using KDE + mutual information

clear all 

samples = 200;

%compile sample 1
mean1 = 3; %mean
set1 = randn(samples,1) + mean1;

%compile sample 2
mean2 = 6;%mean
set2 = randn(samples,1) + mean2;

%Get distributions
bandwidth = 4;
[ prob1, x1 ] = KDE1D( set1, bandwidth );
[ prob2, x2 ] = KDE1D( set2, bandwidth ); %%%%%%%%Why are these different size!?

% compute entropy
entropy1 = EntropyCalc( prob1 );
entropy2 = EntropyCalc( prob2 );

% compute mutual info
%mutinf = mutualInformation(prob1, prob2)

%plot samples
figure(1)
scatter(set1,ones(samples,1)*0.01,'r')
hold on
plot(x1,prob1,'r-','LineWidth',2);
hold on

scatter(set2,ones(samples,1)*0.01,'g')
hold on
plot(x2,prob2,'g-','LineWidth',2);
hold off

xlabel('X (units)')
ylabel('Probability (units)')
title('Test Distributions')
legend('data1','pdf 1','data 2','pdf 2')

%Entropy value plot
figure(2)
plot(1:10,ones(10,1)*entropy1,'r');
hold on
plot(1:10,ones(10,1)*entropy2,'g');
hold off
axis([0,10,0,1.3*max(entropy1,entropy2)])
ylabel('Entropy')
title('Entropy Values')
legend('sample1','sample2')

