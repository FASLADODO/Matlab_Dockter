function Lambdas = Lambda_Critical(data)
%Lambda_Critical(data): Find critical Lambda matrix using tims idea

%Get Sizes
[nn,order] = size(data{1}.state);
classes = length(data);

Mod_States = ModifyStates(data);


Lambdas = [];
%Compute Parameters for each class
for ii = 1:classes
    for jj = 1:classes
        %subtract other classes with lambda
        if ii ~= jj

            %compute lambdars
            Lambdas{ii,jj} = Mod_States{ii}.X_Bar *pinv( Mod_States{jj}.X_Bar  );
        end
    end
end


end