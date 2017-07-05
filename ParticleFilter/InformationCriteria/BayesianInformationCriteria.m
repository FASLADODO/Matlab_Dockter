function [Errors, MeanErrors, Params_All] = BayesianInformationCriteria(X,Y,states,class,option)
%BayesianInformationCritera(X,Y,option): Leave one column of X out each time
%Then compute parameters according to option
%Using parameters compute error relative to Y
%Plot average error for each parameter being left out.
%option= 'TLS'

states = ['AllStates',states];

[NN,SS] = size(X);

Errors = [];
Params_All = [];

for ii = 1:SS
   %Remove one column from data
   X_temp = X;
   if(ii == 1)
       X_temp(:,2) = [];
   end
   X_temp(:,ii) = [];
   params = [];
   
   %compute parameters to use
   if(option == 'TLS')
      params = inv(X_temp'*X_temp)*X_temp'*Y;
   else
      params = inv(X_temp'*X_temp)*X_temp'*Y;
   end
   
   %parameters
   %Params_All(:,ii) = params;
   
   %get errors
   Errors(:,ii) = Y - X_temp*params;
   
end

params = inv(X'*X)*X'*Y;
Error_None = Y - X*params;
Errors = [Error_None, Errors];

%get means
MeanErrors = rms(Errors,1);

%plot that shit
figure
h = bar(MeanErrors);
set(gca,'XTickLabel',states)
xlabel('Left Out State')
ylabel('Mean Error')
sttit = sprintf('Error vs Removed State: %s',class);
title(sttit)



end