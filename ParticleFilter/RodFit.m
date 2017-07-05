function [ estimate ] = RodFit( Class )
%implementing started covariance fit

%class = .data, .classifier;
data = Class.data;
classifier = Class.classifier;

M_nk = zeros(length(unique(classifier)),length(data) );
for ii = 1: length(data)
    M_nk(classifier(ii),ii) = 1;
end

for jj = unique(classifier)
    u(jj) = sum(M_nk(jj,:).*data) / sum(M_nk(jj,:) );
end

