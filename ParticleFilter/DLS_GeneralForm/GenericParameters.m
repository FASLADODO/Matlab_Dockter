function Params = GenericParameters(data)
%Lambda_Critical(data): Find critical Lambda matrix using tims idea

%Get Sizes
[nn,order] = size(data{1}.state);
classes = length(data);

Mod_States = ModifyStates(data);

for ii = 1:classes
    Params(:,ii) = inv(Mod_States{ii}.X_Bar ) * Mod_States{ii}.U_Bar;
end

end