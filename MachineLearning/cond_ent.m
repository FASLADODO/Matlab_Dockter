function result = cond_ent(Y, X)
%https://alliance.seas.upenn.edu/~cis520/wiki/index.php?n=Lectures.DecisionTrees
% Their way
% Calculates the conditional entropy of y given x
result = 0;

tab = tabulate(X);
% Remove zero-entries
tab = tab(tab(:,3)~=0,:);

for i = 1:size(tab,1)
    % Get entropy for y values where x is the current value
    H = EntropySample(Y(X == tab(i,1)));
    % Get probability
    prob = tab(i, 3) / 100;
    % Combine
    result = result + prob * H;
end
