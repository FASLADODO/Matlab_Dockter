function Lambdas = Lambda_Analytical(data)
%Lambda_Analytical(data): Find Lambda matrix using tims idea of X1 X2^{-1}
%Returns pairwise lambda combos

%Get Sizes
[nn,order] = size(data{1}.state);
classes = length(data);

%Compute numeric version of xbar and ubar:
for ii = 1:classes
    for row = 1:order
        for col = 1:order 
            %Same form as symbolic, just easier than subsituting
            X_Bar_num(row,col,ii) = sum( data{ii}.state(:,row).*data{ii}.state(:,col));
        end
    end

end


LambdasAll = [];
%Compute Parameters for each class
for ii = 1:classes
    for jj = 1:classes
        %subtract other classes with lambda
        if ii ~= jj

            %compute lambdars
            %(lambda_{i j})
            LambdasAll{ii,jj} = X_Bar_num(:,:,ii)*pinv(X_Bar_num(:,:,jj));
        end
    end
end

Lambdas = [];
%only take minimum lambdas
for ii = 1:classes
    for jj = 1:classes
        %subtract other classes with lambda
        if ii ~= jj
            %now max
            if( norm(LambdasAll{ii,jj}) > norm(LambdasAll{jj,ii}) )
                
                %change this
                Lambdas{ii,jj} = eye(order) - LambdasAll{jj,ii};
            else
                Lambdas{ii,jj} = LambdasAll{ii,jj};
            end
        end
    end
end

end