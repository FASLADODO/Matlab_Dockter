function Parameters = DLS_TrainGeneral_LambdaMatrix(data, lambda, verbose)
%Computes DLS Parameters for DDA using training data set forarbitrary
%number of classes and states

%Arguments:
%data is a struct containing a substruct for each training class
%Then each class should contain a input vector and a state matrix
%ie:
%data{1}.input = n x 1
%data{1}.state = n x k (k=# states)
%lambda is a vector containing the DLS lambda parameter for each class
%ie:
%lambda = 1 x cc  (cc = # classes)

%Get Sizes
[nn,order] = size(data{1}.state);
classes = length(data);

if(nargin == 2)
    verbose = 0;
end

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
    if(verbose == 1)
        fprintf('class: %i \n',ii)
        fprintf('XMat: ')
        X_Bar_num(:,:,ii)
        fprintf('UMat: ')
        U_Bar_num(:,ii)
    end
end


Parameters = [];
%Compute Parameters for each class using pairwise
for ii = 1:classes
    
    Parameters_0(:,ii) = pinv(X_Bar_num(:,:,ii))*U_Bar_num(:,ii);

    for jj = 1:classes
        %subtract other classes with lambda
        if ii ~= jj
            %should maybe be elementwise
            
            X_mat_sum = X_Bar_num(:,:,ii) - lambda{ii,jj}*X_Bar_num(:,:,jj); %X_1 - L_12*X_2
            U_vec_sum = U_Bar_num(:,ii) -  lambda{ii,jj} *U_Bar_num(:,jj);
            Parameters{ii,jj} = pinv(X_mat_sum)*U_vec_sum;
            
            
            
            if(verbose == 1)
                fprintf('class: %i, opposite: %i \n',ii,jj)
                fprintf('x matrix \n')
                X_mat_sum
                fprintf('U vector \n')
                U_vec_sum
                fprintf('Parameters\n')
                Parameters{ii,jj}
            end
        end
    end
end

Parameters_0

% Parameters


end
