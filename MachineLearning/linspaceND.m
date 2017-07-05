function X = linspaceND(Data,nn)
%creates linearly spaced data in n dimensions

% data is the data vector you wish to span
% Each row is a sample, each column is a dimension

% nn is the number of points to create'


[~,SD] = size(Data);

X = zeros(nn,SD);
for ii = 1:SD
    X(:,ii) = linspace(min(Data(:,ii)),max(Data(:,ii)),nn)';
end