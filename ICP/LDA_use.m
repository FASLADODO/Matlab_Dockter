%http://matlabdatamining.blogspot.com/2010/12/linear-discriminant-analysis-lda.html

%first import data sets into matrix

ppl_per_md = inputdata(:,1);
life_expect = inputdata(:,2);
over_60 = life_expect > 60; %over 60 life expectancy?

%using LDA function
W = LDA(ppl_per_md,over_60);


%calculate linear scores
L = [ones(length(ppl_per_md),1) ppl_per_md] * W';


%compute class probabilities
P = exp(L) ./ repmat(sum(exp(L),2),[1 2]);

plot(P)