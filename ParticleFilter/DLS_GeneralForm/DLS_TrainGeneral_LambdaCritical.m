function Parameters = DLS_TrainGeneral_LambdaCritical(data, lambda, percent)
%Computes DLS Parameters for DDA using training data set and lambda
%critical with percent of ie percent = 0.95

%Get Sizes
[nn,order] = size(data{1}.state);
classes = length(data);

%Compute numeric version of xbar and ubar:
Mod_States = ModifyStates(data);

Parameters = [];
%Compute Parameters for each class using pairwise
for ii = 1:classes

    for jj = 1:classes
        
        %subtract other classes with lambda
        if ii ~= jj
            %should maybe be elementwise
            
            X_mat_sum = Mod_States{ii}.X_Bar -  percent * lambda{ii,jj} * Mod_States{jj}.X_Bar; %X_1 - L_12*X_2
            U_vec_sum = Mod_States{ii}.U_Bar -  percent * lambda{ii,jj} * Mod_States{jj}.U_Bar;
            Parameters{ii,jj} = pinv(X_mat_sum)*U_vec_sum;
            
        end
    end
end


end
