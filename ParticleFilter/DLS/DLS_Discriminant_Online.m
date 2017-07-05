function [ delta1, delta2, class_arr, class_est, class_est_2, noises, alpha ] = DLS_Discriminant_Online( data_in, class_in, ec_1, ec_2, param1_d, param2_d, tindex )
%The discriminant least squares online classification
%data is a struct contaning:
% data{k}.state{1} = x
% data{k}.state{2} = xdot
% data{k}.state{3}= xdotdot
% data.input = system input
% ie column vectors of all state data
% returns -1 if no info

if(tindex == 0)
    tindex = length(data_in{class_in}.state{1});
end

%online data matrix
Dx = [data_in{class_in}.state{1},data_in{class_in}.state{2},data_in{class_in}.state{3}];

%all error values DIS IS THE SHIT
error_1 = abs( data_in{class_in}.input - (Dx*param1_d) );
error_2 = abs( data_in{class_in}.input - (Dx*param2_d) );

class = 0;
class_est = 0;
conf_sum_1 = [];
conf_sum_2 = [];
cnt_a = 1;
cnt_1 = 1;
cnt_2 = 1;
threshm = 2;
alpha = [];

for tt = 11:tindex

    %up till time t error sum
    delta1 = sum( error_1(tt-10:tt) );
    delta2 = sum( error_2(tt-10:tt) );
    
    error_sum_class_1 = sum(ec_1(tt-10:tt));
    error_sum_class_2 = sum(ec_2(tt-10:tt));
    
    %confidence value
    alpha(tt,1) = (threshm*abs(delta2-delta1)) / (error_sum_class_1 + error_sum_class_2);
    %alpha(tt,1) = (threshm*abs(error_2(tt)-error_1(tt))) / (ec_1(tt) + ec_2(tt));
    
    %compare deltas
    if(delta1 < delta2)
        class_arr(cnt_a) = 1;
        conf_sum_1(cnt_1) = alpha(tt,1);
        cnt_1 = cnt_1 + 1;
    elseif(delta2 < delta1)
        class_arr(cnt_a) = 2;
        conf_sum_2(cnt_2) = alpha(tt,1);
        cnt_2 = cnt_2 + 1;
    else
        class_arr(cnt_a) = 0;
    end
    
    noises(cnt_a,:) = [error_1(tt), error_2(tt), ec_1(tt), ec_2(tt), delta1, delta2, error_sum_class_1, error_sum_class_2 ];

    cnt_a = cnt_a + 1;
end

class_est = mode(class_arr);
class_est_2 = 0;
if(sum(conf_sum_1) > sum(conf_sum_2) )
    class_est_2 = 1;
elseif (sum(conf_sum_2) > sum(conf_sum_1) )
    class_est_2 = 2;
else
    class_est_2 = class_est;
end

end

