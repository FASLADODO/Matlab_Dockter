


%Get Sizes
order = 3;
classes = 2;
nn = 10;

lambda = [0.5,0.5];

data{1}.state = [ones(nn,1)*1,ones(nn,1)*2,ones(nn,1)*3];
data{1}.input = ones(nn,1)*4;

data{2}.state = [ones(nn,1)*5,ones(nn,1)*6,ones(nn,1)*7];
data{2}.input = ones(nn,1)*9;


%%%%%%%%%%%%% SYMBOLIC STUFF %%%%%%%%%%%%%%%%%%%

%create array of symbolic X's
XS = sym(zeros(1, order));
for k=1:order
    XS(k) = sym(sprintf('x%d', k));
end


%Create symbolic U vector
syms U
U_BAR = sym(zeros(order, 1));
for row = 1:order
    U_BAR(row,1) = XS(row)*U;
end

%Create symbolic lambdas 
lambda_sym = sym(zeros(1, classes));
for ii = 1:classes
    lambda_sym(ii) = sym(sprintf('L%d', ii));
end

%Create generic XBar data matrix, This will become a sum over time
%See Equation 2.0, Page 111 in Lab Notebook
X_BAR = sym(zeros(order, order));
for row = 1:order
    for col = 1:order 
        X_BAR(row,col) = XS(row)*XS(col);
    end
end

%Compute Form

%see equation 1.9, page 110 Rods Lab Notebook
X_sym_sum = X_BAR;
U_sym_sum = U_BAR;
for jj = 2:classes
    %subtract other classes with lambda
    X_sym_sum = X_sym_sum - lambda_sym(jj)*X_BAR;
    U_sym_sum = U_sym_sum - lambda_sym(jj)*U_BAR;
end
% X_sym_sum
% U_sym_sum



%%%%%%%%%%%%%%%%%%%% NUMERIC STUFF %%%%%%%%%%%%%%%%%%%%%

%Compute numeric version of xbar and ubar:
for ii = 1:classes
    for row = 1:order
        for col = 1:order 
            %Same form as symbolic, just easier than subsituting
            X_Bar_num(row,col,ii) = sum( data{ii}.state(:,row).*data{ii}.state(:,col));
        end
    end
    
    for row = 1:order
        %Same form as symbolic, just easier than subsituting
        U_Bar_num(row,ii) = sum( data{ii}.state(:,row).*data{ii}.input );
    end

end


Parameters = [];
%Compute Parameters for each class
for ii = 1:classes
    %see equation 1.9, page 110 Rods Lab Notebook
    X_mat_sum = X_Bar_num(:,:,ii);
    U_vec_sum = U_Bar_num(:,ii);
    for jj = 1:classes
        %subtract other classes with lambda
        if ii ~= jj
            X_mat_sum = X_mat_sum - lambda(jj)*X_Bar_num(:,:,jj);
            U_vec_sum = U_vec_sum - lambda(jj)*U_Bar_num(:,jj);
        end
    end
    X_mat_sum
    U_vec_sum
    %Compute DLS Parameters for each class
    Parameters(:,ii) = pinv(X_mat_sum)*U_vec_sum;
end

Parameters

%% Other Way

data{1}.state = [];
data{2}.state = [];

data{1}.state{1} = ones(nn,1)*1;
data{1}.state{2} = ones(nn,1)*2;
data{1}.state{3} = ones(nn,1)*3;

data{2}.state{1} = ones(nn,1)*5;
data{2}.state{2} = ones(nn,1)*6;
data{2}.state{3} = ones(nn,1)*7;




lambda = 0.5;


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
        m31,m32,m33]
vec1 = [v1;v2;v3]
param1 = pinv(mat1) * vec1


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
        m31,m32,m33]
vec2 = [v1;v2;v3]
param2 = pinv(mat2) * vec2




