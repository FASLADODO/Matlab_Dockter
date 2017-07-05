function [ param1, param2 ] = DLS_Discriminant_Train( data, lambda )
%The discriminant least squares (for binary classification) only 3rd order
%systems
%data is a struct contaning:
% data{k}.state{1} = x
% data{k}.state{2} = xdot
% data{k}.state{3}= xdotdot
% data.input = system input
% ie column vectors of all state data
% Where {k} indicates the corresponding originating class label
%lambda is a weighting factor between generative and discirminative models

%Returns: paramvectors 1,2,d = [K;D;M]


%This assumes inputs are same length


%compute all time summations for inverse matrix


m11 = sum( data{1}.state{1}.^2 - lambda*(data{2}.state{1}.^2) );
m12 = sum( data{1}.state{1}.*data{1}.state{2} - lambda*(data{2}.state{1}.*data{2}.state{2}) );
m13 = sum( data{1}.state{1}.*data{1}.state{3} - lambda*(data{2}.state{1}.*data{2}.state{3}) );

m21 = sum( data{1}.state{1}.*data{1}.state{2} - lambda*(data{2}.state{1}.*data{2}.state{2}) ); 
m22 = sum( data{1}.state{2}.^2 - lambda*(data{2}.state{2}.^2) );
m23 = sum( data{1}.state{2}.*data{1}.state{3} - lambda*(data{2}.state{2}.*data{2}.state{3}) );

m31 = sum( data{1}.state{1}.*data{1}.state{3} - lambda*(data{2}.state{1}.*data{2}.state{3}) ); 
m32 = sum( data{1}.state{2}.*data{1}.state{3} - lambda*(data{2}.state{2}.*data{2}.state{3}) ); 
m33 = sum( data{1}.state{3}.^2 - lambda*(data{2}.state{3}.^2) );

%compute time summations with inputs
v1 = sum(data{1}.state{1}.*data{1}.input - (lambda*data{2}.state{1}.*data{2}.input) );
v2 = sum(data{1}.state{2}.*data{1}.input - (lambda*data{2}.state{2}.*data{2}.input) );
v3 = sum(data{1}.state{3}.*data{1}.input - (lambda*data{2}.state{3}.*data{2}.input) );


%Compute parameters 1 vector
mat1 = [m11,m12,m13;
        m21,m22,m23;
        m31,m32,m33];
vec1 = [v1;v2;v3];
param1 = pinv(mat1) * vec1;


%now compute second params
% all time summations for inverse matrix

m11 = sum( data{2}.state{1}.^2 - lambda*(data{1}.state{1}.^2) );
m12 = sum( data{2}.state{1}.*data{2}.state{2} - lambda*(data{1}.state{1}.*data{1}.state{2}) );
m13 = sum( data{2}.state{1}.*data{2}.state{3} - lambda*(data{1}.state{1}.*data{1}.state{3}) );

m21 = sum( data{2}.state{1}.*data{2}.state{2} - lambda*(data{1}.state{1}.*data{1}.state{2}) ); 
m22 = sum( data{2}.state{2}.^2 - lambda*(data{1}.state{2}.^2) );
m23 = sum( data{2}.state{2}.*data{2}.state{3} - lambda*(data{1}.state{2}.*data{1}.state{3}) );

m31 = sum( data{2}.state{1}.*data{2}.state{3} - lambda*(data{1}.state{1}.*data{1}.state{3}) ); 
m32 = sum( data{2}.state{2}.*data{2}.state{3} - lambda*(data{1}.state{2}.*data{1}.state{3}) ); 
m33 = sum( data{2}.state{3}.^2 - lambda*(data{1}.state{3}.^2) );

%compute time summations with inputs
v1 = sum(data{2}.state{1}.*data{2}.input - (lambda*data{1}.state{1}.*data{1}.input) );
v2 = sum(data{2}.state{2}.*data{2}.input - (lambda*data{1}.state{2}.*data{1}.input) );
v3 = sum(data{2}.state{3}.*data{2}.input - (lambda*data{1}.state{3}.*data{1}.input) );

%Compute parameters 1 vector
mat2 = [m11,m12,m13;
        m21,m22,m23;
        m31,m32,m33];
vec2 = [v1;v2;v3];
param2 = pinv(mat2) * vec2;




end

