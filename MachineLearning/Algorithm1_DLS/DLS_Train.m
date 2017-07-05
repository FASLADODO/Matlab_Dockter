function Parameters = DLS_Train(X, Y, Labels, lambda)
%Computes DLS Parameters for DDA using training data set forarbitrary
%number of classes and states

%Arguments:
%X is the data matrix for all classes
%Y is the output matrix for all classes
%Labels is the class labels
%lambda is a vector containing the DLS lambda parameter for each class
%ie:
%lambda = 1 x cc  (cc = # classes)

%Get Sizes
[NN,SS] = size(X);
cslist = unique(Labels);
classes = length(cslist);

%make sure we guci
if(length(lambda) ~= classes)
    lambda = repmat(lambda(1),1,classes);
end

%big ol storage structs
data = [];

%get sub sample data for each class
for cc = 1:classes
    tempX = X(Labels == cslist(cc),:);
    tempY = Y(Labels == cslist(cc),:);
    
    %store into struct
    data{cc}.state = tempX;
    data{cc}.input = tempY;
end

%Compute numeric version of xbar and ubar:
for cc = 1:classes
    X_Bar_num(:,:,cc) = data{cc}.state' *data{cc}.state;
    U_Bar_num(:,cc) = data{cc}.state'*data{cc}.input;
end      

Parameters = [];
%Compute Parameters for each class
for cc = 1:classes
    %see equation 1.9, page 110 Rods Lab Notebook
    X_mat_sum = X_Bar_num(:,:,cc);
    U_vec_sum = U_Bar_num(:,cc);
    for jj = 1:classes
        %subtract other classes with lambda
        if cc ~= jj
            X_mat_sum = X_mat_sum - lambda(jj)*X_Bar_num(:,:,jj);
            U_vec_sum = U_vec_sum - lambda(jj)*U_Bar_num(:,jj);
        end
    end
    %Compute DLS Parameters for each class
    Parameters(:,cc) = pinv(X_mat_sum)*U_vec_sum;
end


end
