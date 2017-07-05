function [P,Class] = LDAonline(X,W,C)
%using LDA weights to compute class probabilites
% Use:
% [W,C] = LDA(X,Y,Priors)
%
% P     = Class Probabilities, eacolumns are classes, rows are samples
% Class     = greatest probabilities class estimate
% X   = predictor data (variables in columns, observations in rows)
% W  = LDA coefficients, each row is for a class
% C     = constants, each row is for a class


% observation and weights size
[NN,~] = size(X);
[NW,~] = size(W);

% Calulcate linear scores for training data
L =  ( W*X' + repmat(C,1,NN) )';

% Calculate class probabilities
P = exp(L) ./ repmat( sum( exp(L),2 ), [1 NW] );

%find class estimate for each sample
[~,Class] = max(P,[],2);


end