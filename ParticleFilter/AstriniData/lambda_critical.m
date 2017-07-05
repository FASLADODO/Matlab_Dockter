function lambdas = lambda_critical(data)
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


lambdas = [];
%Compute Parameters for each class
for ii = 1:classes

    for jj = 1:classes
        %subtract other classes with lambda and scaled by difference in
        %data points
        if ii ~= jj
            lambda.class{ii}.cross{jj} =  X_Bar_num(:,:,ii)*inv(X_Bar_num(:,:,jj));
            %lambdas(ii,jj) =  det( X_Bar_num(:,:,ii)*inv(X_Bar_num(:,:,jj)) );
        end
    end
end


end