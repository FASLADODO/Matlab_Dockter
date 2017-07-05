function [ gamma ] = Compute_Gamma( parambar, param )
%Compute gamma, parameter similarity measure

gamma = norm(parambar - param)/norm(param);

end
