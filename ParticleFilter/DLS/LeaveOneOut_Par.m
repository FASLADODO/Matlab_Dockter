function [ classification ] = LeaveOneOut_Par( pp, simd, lambdas, runs, class_actual, ec_1, ec_2 )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    for L = 1:length(lambdas)
        
        for rr = 1:runs
            cnt_all=(pp-1)*runs+rr;
            %setup indexes
            idx_train = 1:runs;
            idx_train(idx_train == rr) = [];
            
            %Clear?
            simd{1}.state{1} = [];
            simd{2}.state{2} = [];
            simd{1}.state{3} = [];
            simd{1}.input = [];
            simd{2}.state{1} = [];
            simd{2}.state{2} = [];
            simd{2}.state{3} = [];
            simd{2}.input = [];
            
            simd{1}.state{1} = reshape(simd{1}.state_k{1}(:,idx_train),[],1);
            simd{1}.state{2} = reshape(simd{1}.state_k{2}(:,idx_train),[],1);
            simd{1}.state{3} = reshape(simd{1}.state_k{3}(:,idx_train),[],1);
            simd{1}.input = reshape(simd{1}.input_k(:,idx_train),[],1);
            
            simd{2}.state{1} = reshape(simd{2}.gamma{pp}.state_k{1}(:,idx_train),[],1);
            simd{2}.state{2} = reshape(simd{2}.gamma{pp}.state_k{2}(:,idx_train),[],1);
            simd{2}.state{3} = reshape(simd{2}.gamma{pp}.state_k{3}(:,idx_train),[],1);
            simd{2}.input = reshape(simd{2}.gamma{pp}.input_k(:,idx_train),[],1);
            
            %resample online data
            if(class_actual == 2)
                online{class_actual}.state{1} = reshape(simd{class_actual}.gamma{pp}.state_k{1}(:,rr),[],1);
                online{class_actual}.state{2} = reshape(simd{class_actual}.gamma{pp}.state_k{2}(:,rr),[],1);
                online{class_actual}.state{3} = reshape(simd{class_actual}.gamma{pp}.state_k{3}(:,rr),[],1);
                online{class_actual}.input = reshape(simd{class_actual}.gamma{pp}.input_k(:,rr),[],1);
            else
                online{class_actual}.state{1} = reshape(simd{class_actual}.state_k{1}(:,rr),[],1);
                online{class_actual}.state{2} = reshape(simd{class_actual}.state_k{2}(:,rr),[],1);
                online{class_actual}.state{3} = reshape(simd{class_actual}.state_k{3}(:,rr),[],1);
                online{class_actual}.input = reshape(simd{class_actual}.input_k(:,rr),[],1);
            end
            
            
            %str_info = sprintf('TRAIN L =  %f, G = %f, R = %f \n',lambdas(L),gamma(pp),rr);
            %Perform training
            %fprintf(str_info)
            [ param1(:,cnt_all), param2(:,cnt_all) ] = DLS_Discriminant_Train( simd, lambdas(L) );
            
            dummy{1}.state = [simd{1}.state{3}, simd{1}.state{2}, simd{1}.state{1}];
            dummy{2}.state = [simd{2}.state{3}, simd{2}.state{2}, simd{2}.state{1}];
            dummy{1}.input = simd{1}.input;
            dummy{2}.input = simd{2}.input;
            
            GParams = DLS_TrainGeneral(dummy, [lambdas(L),lambdas(L)] );
            
            param1_G(:,cnt_all) = flipud(GParams(:,1));
            param2_G(:,cnt_all) = flipud(GParams(:,2));
            
            %do online classification
            %fprintf(' ONLINE\n')
            [ delta1, delta2, class_arr, class_est, class_est_2, noises_vec, alpha ] = DLS_Discriminant_Online( online, class_actual,  ec_1, ec_2(:,pp), param1(:,cnt_all), param2(:,cnt_all), 0 );
            
            %store data
            classification.lambda{L}.class(:,rr) = [class_est;class_actual;class_est_2];
            classification.lambda{L}.nses{rr} = noises_vec;
            classification.lambda{L}.class_arr(:,rr) = class_arr';
            classification.lambda{L}.alpha(:,rr)= alpha;
            classification.lambda{L}.params(:,rr) = [param1(:,cnt_all);param2(:,cnt_all)];
        end
    end

end

