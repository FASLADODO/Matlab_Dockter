function [ params ] = TotalLeastSquares( data, class )
%Normal total least squares

%data is a struct contaning:
% data{k}.state{1} = x
% data{k}.state{2} = xdot
% data{k}.state{3}= xdotdot
% data.input = system input
% ie column vectors of all state data
% K represents class

if(length(data{class}.state{1}) ~= length(data{class}.input) )
    fprintf('states and inputs are different lengths') 
end

%make D matrix
D = [data{class}.state{1},data{class}.state{2},data{class}.state{3} ];

%Tims trick to estimate phi
params = inv(D' * D) * D' * data{class}.input;

end

