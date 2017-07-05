%test 1d data idea

%TODO, use tims notes

nn = 1000;
fontS = 13;

center_1 = [5,7];
center_2 = [5,6.5];
spread = 0.5; %2

set1 = [randn(nn,1).*spread + center_1(1), randn(nn,1).*spread + center_1(2)]; %around 5,0
set2 = [randn(nn,1).*spread + center_2(1), randn(nn,1).*spread + center_2(2)]; %around 3,-3

figure
h1 = plot(set1(:,1),set1(:,2),'r.','MarkerSize',5);
hold on
h2 = plot(set2(:,1),set2(:,2),'g.','MarkerSize',5);

hold off
title('Distributions for Two Classes','FontSize',fontS)
xlabel('U','FontSize',fontS)
ylabel('X','FontSize',fontS)


%in DLS form
X = [set1(:,2).*set1(:,2),set2(:,2).*set2(:,2)];
U = [set1(:,1).*set1(:,2),set2(:,1).*set2(:,2)];

%get lambda critical
lambda_critical(1) = sum(X(:,1))*inv(sum(X(:,2)));
lambda_critical(2) = sum(X(:,2))*inv(sum(X(:,1)));

%%

lambda_percent = 0:0.01:0.99;


%get params
for lp = 1:length(lambda_percent)
    params1(lp) = inv( sum(X(:,1)) - lambda_percent(lp)*lambda_critical(1)*sum(X(:,2)) ) * ( sum(U(:,1)) - lambda_percent(lp)*lambda_critical(1)*sum(U(:,2)) ); 
    params2(lp) = inv( sum(X(:,2)) - lambda_percent(lp)*lambda_critical(2)*sum(X(:,1)) ) * ( sum(U(:,2)) - lambda_percent(lp)*lambda_critical(2)*sum(U(:,1)) ); 
    
end

%general params
param1_g = inv(sum(X(:,1))) *sum(U(:,1))
param2_g = inv(sum(X(:,2))) *sum(U(:,2))

param1_kn = center_1(1) / center_1(2)
param2_kn = center_2(1) / center_2(2)


%plot params 1
figure
scatter( lambda_percent, params1, 'r.');
hold on
plot( [lambda_percent(1),lambda_percent(end)], [param1_g,param1_g], 'b' );

hold off
title('Change in Parameter Class 1','FontSize',fontS)
xlabel('\zeta','FontSize',fontS)
ylabel('Phi','FontSize',fontS)
legend('DLS params','true parameter')

%plot params 2
figure
scatter( lambda_percent, params2, 'r.');
hold on
plot( [lambda_percent(1),lambda_percent(end)], [param2_g,param2_g], 'b' );

hold off
title('Change in Parameter Class 2','FontSize',fontS)
xlabel('\zeta','FontSize',fontS)
ylabel('Phi','FontSize',fontS)
legend('DLS param','true parameter')




%get params
for ii = 1:2
    for lp = 1:length(lambda_percent)
        %get errors at each time step
        error1.classActual{ii}(:,lp) = abs( U(:,ii) - params1(lp)*X(:,ii) );
        error2.classActual{ii}(:,lp) = abs( U(:,ii) - params2(lp)*X(:,ii) );
        
        %get cumsum errors
        sumerror.classActual{ii}.zeta{lp} = [cumsum(error1.classActual{ii}(:,lp)), cumsum(error2.classActual{ii}(:,lp))];
        
        %gettin min
        [val,idx] = min( sumerror.classActual{ii}.zeta{lp}, [],2);
        
        %make class version
        classify.actual{ii}.zeta{lp} = [idx, ones(length(idx),1)*ii ];
        
        corrects = classify.actual{ii}.zeta{lp}(:,1) == classify.actual{ii}.zeta{lp}(:,2);
        percentage(ii,lp) = sum(corrects)/length(corrects);
    end
end


%%

% plot accuracy
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

