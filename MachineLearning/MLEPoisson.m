%%MLE Poisson

x = 0:15;
lambda = 5
noiz = 0.001;

%get probabilities
fx = (exp(-lambda).*lambda.^(x))./(factorial(x)) + randn(1,length(x))*noiz;

figure
plot(x,fx,'+')


%Now get discrete versions
nn = 1000;
sample = [];
for id = 1:length(x)
    sample = [sample, ones(1,round(fx(id)*nn))*x(id)]; %x = 1;
    vals(id) = round(fx(id)*nn);
end

figure
plot(x,vals,'d')

%Maximize Ln(L)
%http://www.colorado.edu/economics/morey/7818/estimation/maxlik/maxlik.pdf

xbar = mean(sample)

L = 0:0.01:10;

Lstar = -L + log(L).*xbar;

[m,idx] = max(Lstar);

Lml = L(idx)

figure
plot(L,Lstar)
hold on
scatter(L(idx),m,'*')
hold off

