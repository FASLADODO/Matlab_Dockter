function Parameters = DLS_TrainGeneral(data, lambda, symbolic)
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
   symbolic = 0; 
end

if symbolic == 1
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
    X_sym_sum
    U_sym_sum

end


%%%%%%%%%%%%%%%%%%%% NUMERIC STUFF %%%%%%%%%%%%%%%%%%%%%

%Get number of data points in each class
for ii = 1:classes
    n_dp(ii) = length(data{ii}.state);
    n_dp_input(ii) = length(data{ii}.input);
    if(n_dp(ii) ~= n_dp_input(ii))
       error('state and input must be same length!') 
    end
end

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
        %subtract other classes with lambda and scaled by difference in
        %data points
        if ii ~= jj
            X_mat_sum = X_mat_sum - (lambda(jj) * (n_dp(ii)/n_dp(jj)) * X_Bar_num(:,:,jj));
            U_vec_sum = U_vec_sum - (lambda(jj) * (n_dp(ii)/n_dp(jj)) * U_Bar_num(:,jj));
%             X_mat_sum = X_mat_sum - lambda(jj) *  X_Bar_num(:,:,jj);
%             U_vec_sum = U_vec_sum - lambda(jj) *  U_Bar_num(:,jj);
        end
    end
%     X_mat_sum
%     U_vec_sum
    %Compute DLS Parameters for each class
    Parameters(:,ii) = pinv(X_mat_sum)*U_vec_sum;
end

% Parameters


end
