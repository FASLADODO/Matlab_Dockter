%test 1d data idea

%TODO, use tims notes

nn = 500;
fontS = 13;

center_1 = [5,7];
center_2 = [5,9];
center_3 = [7,9];
spread = 1; %2

set1 = [randn(nn,1).*spread + center_1(1), randn(nn,1).*spread + center_1(2)]; 
set2 = [randn(nn,1).*spread + center_2(1), randn(nn,1).*spread + center_2(2)]; 
set3 = [randn(nn,1).*spread + center_3(1), randn(nn,1).*spread + center_3(2)];

figure
h1 = plot(set1(:,1),set1(:,2),'r.','MarkerSize',5);
hold on
h2 = plot(set2(:,1),set2(:,2),'g.','MarkerSize',5);
hold on
h3 = plot(set3(:,1),set3(:,2),'b.','MarkerSize',5);
hold off
title('Distributions for Two Classes','FontSize',fontS)
xlabel('U','FontSize',fontS)
ylabel('X','FontSize',fontS)


%in DLS form
X = [set1(:,2).*set1(:,2),set2(:,2).*set2(:,2),set3(:,2).*set3(:,2)];
U = [set1(:,1).*set1(:,2),set2(:,1).*set2(:,2),set3(:,1).*set3(:,2)];

%get lambda critical
lambda_critical(1,2) = sum(X(:,1))*inv(sum(X(:,2)));
lambda_critical(1,3) = sum(X(:,1))*inv(sum(X(:,3)));
lambda_critical(2,1) = sum(X(:,2))*inv(sum(X(:,1)));
lambda_critical(2,3) = sum(X(:,2))*inv(sum(X(:,3)));
lambda_critical(3,1) = sum(X(:,3))*inv(sum(X(:,1)));
lambda_critical(3,2) = sum(X(:,3))*inv(sum(X(:,2)));

disp('lambda critical: ')
lambda_critical

%%

lambda_percent = 0:0.01:0.99;


%get params
for lp = 1:length(lambda_percent)
    %lambda percent
    params{1}{2}(lp) = inv( sum(X(:,1)) - lambda_percent(lp)*lambda_critical(1,2)*sum(X(:,2)) ) * ( sum(U(:,1)) - lambda_percent(lp)*lambda_critical(1,2)*sum(U(:,2)) ); 
    params{1}{3}(lp) = inv( sum(X(:,1)) - lambda_percent(lp)*lambda_critical(1,3)*sum(X(:,3)) ) * ( sum(U(:,1)) - lambda_percent(lp)*lambda_critical(1,3)*sum(U(:,3)) ); 
    params{2}{1}(lp) = inv( sum(X(:,2)) - lambda_percent(lp)*lambda_critical(2,1)*sum(X(:,1)) ) * ( sum(U(:,2)) - lambda_percent(lp)*lambda_critical(2,1)*sum(U(:,1)) ); 
    params{2}{3}(lp) = inv( sum(X(:,2)) - lambda_percent(lp)*lambda_critical(2,3)*sum(X(:,3)) ) * ( sum(U(:,2)) - lambda_percent(lp)*lambda_critical(2,3)*sum(U(:,3)) ); 
    params{3}{1}(lp) = inv( sum(X(:,3)) - lambda_percent(lp)*lambda_critical(3,1)*sum(X(:,1)) ) * ( sum(U(:,3)) - lambda_percent(lp)*lambda_critical(3,1)*sum(U(:,1)) ); 
    params{3}{2}(lp) = inv( sum(X(:,3)) - lambda_percent(lp)*lambda_critical(3,2)*sum(X(:,2)) ) * ( sum(U(:,3)) - lambda_percent(lp)*lambda_critical(3,2)*sum(U(:,2)) ); 
end

%general params
param1_g = inv(sum(X(:,1))) *sum(U(:,1))
param2_g = inv(sum(X(:,2))) *sum(U(:,2))
param3_g = inv(sum(X(:,3))) *sum(U(:,3))

param1_kn = center_1(1) / center_1(2)
param2_kn = center_2(1) / center_2(2)
param3_kn = center_3(1) / center_3(2)


%plot params 1 - 2
figure
scatter( lambda_percent, params{1}{2}, 'r.');
hold on
plot( [lambda_percent(1),lambda_percent(end)], [param1_g,param1_g], 'b' );

hold off
title('Change in Parameter Class 1-2','FontSize',fontS)
xlabel('\zeta','FontSize',fontS)
ylabel('Phi','FontSize',fontS)
legend('DLS params','true parameter')

%plot params 1 - 3
figure
scatter( lambda_percent, params{1}{3}, 'r.');
hold on
plot( [lambda_percent(1),lambda_percent(end)], [param1_g,param1_g], 'b' );

hold off
title('Change in Parameter Class 1-3','FontSize',fontS)
xlabel('\zeta','FontSize',fontS)
ylabel('Phi','FontSize',fontS)
legend('DLS params','true parameter')

%plot params 2-1
figure
scatter( lambda_percent, params{2}{1}, 'r.');
hold on
plot( [lambda_percent(1),lambda_percent(end)], [param2_g,param2_g], 'b' );

hold off
title('Change in Parameter Class 2-1','FontSize',fontS)
xlabel('\zeta','FontSize',fontS)
ylabel('Phi','FontSize',fontS)
legend('DLS param','true parameter')


%plot params 2-3
figure
scatter( lambda_percent, params{2}{3}, 'r.');
hold on
plot( [lambda_percent(1),lambda_percent(end)], [param2_g,param2_g], 'b' );

hold off
title('Change in Parameter Class 2-3','FontSize',fontS)
xlabel('\zeta','FontSize',fontS)
ylabel('Phi','FontSize',fontS)
legend('DLS param','true parameter')


%plot params 3-1
figure
scatter( lambda_percent, params{3}{1}, 'r.');
hold on
plot( [lambda_percent(1),lambda_percent(end)], [param3_g,param3_g], 'b' );

hold off
title('Change in Parameter Class 3-1','FontSize',fontS)
xlabel('\zeta','FontSize',fontS)
ylabel('Phi','FontSize',fontS)
legend('DLS param','true parameter')

%plot params 3-2
figure
scatter( lambda_percent, params{3}{2}, 'r.');
hold on
plot( [lambda_percent(1),lambda_percent(end)], [param3_g,param3_g], 'b' );

hold off
title('Change in Parameter Class 3-2','FontSize',fontS)
xlabel('\zeta','FontSize',fontS)
ylabel('Phi','FontSize',fontS)
legend('DLS param','true parameter')



%% compute errors

lc = 0.95;

for ii = 1:3
    for lp = 1:length(lambda_percent)
        %get errors at each time step
        for jj = 1:3
            for kk = 1:3
                if( jj ~= kk)
                    error{jj}{kk}.classActual{ii}(:,lp) = abs( U(:,ii) - params{jj}{kk}(lp)*X(:,ii) );
                    
                    sumerror{jj}{kk}.classActual{ii}.zeta{lp} = cumsum(error{jj}{kk}.classActual{ii}(:,lp));
                    
                    if ( lambda_percent(lp) == lc )
                        classing{ii}(jj,kk) = sum(error{jj}{kk}.classActual{ii}(:,lp));
                    end
                end
            end
        end
        
%         temp1 = ( sumerror{1}{2}.classActual{ii}.zeta{lp} + sumerror{1}{3}.classActual{ii}.zeta{lp} )./( sumerror{2}{1}.classActual{ii}.zeta{lp} + sumerror{3}{1}.classActual{ii}.zeta{lp} );
%         temp2 = ( sumerror{2}{1}.classActual{ii}.zeta{lp} + sumerror{2}{3}.classActual{ii}.zeta{lp} )./( sumerror{1}{2}.classActual{ii}.zeta{lp} + sumerror{3}{2}.classActual{ii}.zeta{lp} );
%         temp3 = ( sumerror{3}{1}.classActual{ii}.zeta{lp} + sumerror{3}{2}.classActual{ii}.zeta{lp} )./( sumerror{1}{3}.classActual{ii}.zeta{lp} + sumerror{2}{3}.classActual{ii}.zeta{lp} );
%         
        temp1 = ( sumerror{1}{2}.classActual{ii}.zeta{lp} ./ sumerror{2}{1}.classActual{ii}.zeta{lp} ) + ( sumerror{1}{3}.classActual{ii}.zeta{lp} ./ sumerror{3}{1}.classActual{ii}.zeta{lp} );
        temp2 = ( sumerror{2}{1}.classActual{ii}.zeta{lp} ./ sumerror{1}{2}.classActual{ii}.zeta{lp} ) + ( sumerror{2}{3}.classActual{ii}.zeta{lp} ./ sumerror{3}{2}.classActual{ii}.zeta{lp} );
        temp3 = ( sumerror{3}{1}.classActual{ii}.zeta{lp} ./ sumerror{1}{3}.classActual{ii}.zeta{lp} ) + ( sumerror{3}{2}.classActual{ii}.zeta{lp} ./ sumerror{2}{3}.classActual{ii}.zeta{lp} );
        
        
        ErrorRatio{ii}.zeta{lp} = [temp1,temp2,temp3];
        
        %gettin min
        [val,idx] = min( ErrorRatio{ii}.zeta{lp}, [],2);
        
        %make class version
        classify.actual{ii}.zeta{lp} = [idx, ones(length(idx),1)*ii ];
        
        corrects = classify.actual{ii}.zeta{lp}(:,1) == classify.actual{ii}.zeta{lp}(:,2);
        percentage(ii,lp) = sum(corrects)/length(corrects);
    end
end


%%

% plot accuracies
figure
cls = 1;
plot( lambda_percent,percentage(cls,:), 'b' );

hold off
title('Accuracy Vs % of \lambda Critical (Class 1)','FontSize',fontS)
xlabel('% \lambda Critical','FontSize',fontS)
ylabel('Accuracy','FontSize',fontS)

figure
cls = 2;
plot( lambda_percent,percentage(cls,:), 'b' );

hold off
title('Accuracy Vs % of \lambda Critical (Class 2)','FontSize',fontS)
xlabel('% \lambda Critical','FontSize',fontS)
ylabel('Accuracy','FontSize',fontS)

figure
cls = 3;
plot( lambda_percent,percentage(cls,:), 'b' );

hold off
title('Accuracy Vs % of \lambda Critical (Class 3)','FontSize',fontS)
xlabel('% \lambda Critical','FontSize',fontS)
ylabel('Accuracy','FontSize',fontS)


