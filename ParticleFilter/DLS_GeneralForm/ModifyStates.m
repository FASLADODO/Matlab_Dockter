function Mod_States = ModifyStates(data)
%Lambda_Analytical(data): Get modified states X_Bar and U_Bar
%state and input are in their own structs

%Get Sizes
[nn,order] = size(data{1}.state);
classes = length(data);

%Compute numeric version of xbar and ubar:
for ii = 1:classes
    for row = 1:order
        for col = 1:order 
            %Same form as symbolic, just easier than subsituting
            %Mod_States{ii}.X_Bar(row,col) = sum( data{ii}.state(:,row).*data{ii}.state(:,col));
            Mod_States{ii}.X_Bar(row,col) = sum( data{ii}.state{row}.*data{ii}.state{col});
        end
    end

end

for ii = 1:classes
    for row = 1:order
        %Same form as symbolic, just easier than subsituting
        %Mod_States{ii}.U_Bar(row,1) = sum( data{ii}.state(:,row).*data{ii}.input );
        Mod_States{ii}.U_Bar(row,1) = sum( data{ii}.state{row}.*data{ii}.input );
    end
end



end