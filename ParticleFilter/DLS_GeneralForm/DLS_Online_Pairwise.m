function Classify = DLS_Online_Pairwise(data1, data2, class_online, parameters)
%%% TODOOOO

%REMOVE THIS LATER PROBABLY
data{1} = data1;
data{2} = data2;


%Get Sizes
[nn,order] = size(data{1}.state);
classes = length(data);


D_online = data{class_online}.state;
U_online = data{class_online}.input;

min_ind = length(U_online);


%Loop through all classes, get errors
for ii = 1:classes
    
    %loop through all pairwise things
    for jj = 1:classes
        if ii ~= jj
            delta{ii,jj} = abs( U_online - (D_online * parameters{ii,jj} ) );
            sumcum{ii,jj} = cumsum(delta{ii,jj},1);
        end
    end
end


%loop through time steps
for kk = 1:min_ind

    %Get classes from errors
    for ii = 1:classes
        num_ratio = 0;
        den_ratio = 0;

        %loop through all pairwise things
        for jj = 1:classes
            if ii ~= jj
                num_ratio = num_ratio + sumcum{ii,jj}(kk);
                den_ratio = den_ratio + sumcum{jj,ii}(kk);
            end
        end
        
        %compute ratios
        if(den_ratio ~= 0 )
            Classify.ratio(kk,ii) = num_ratio/den_ratio;
        end
        
        %Save numerators/denominators
        Classify.num(kk,ii) = num_ratio;
        Classify.den(kk,ii) = den_ratio;
    end
    
    [min_ratio,class_est] = min(Classify.ratio(kk,:) );
    
    Classify.classes(kk,:) = [class_est,class_online];
end

end