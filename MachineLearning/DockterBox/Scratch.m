[valnov,inov] = sort(mean(storeNOVmean));
[valexp,iexp] = sort(storeEXPmean,'descend');
storeNOVstd = std(storeNOVmean);

figure
errorbar(1:length(valnov),100*valnov,100*storeNOVstd(inov))
hold on
% figure
% plot(1:length(valexp),100*valexp)
errorbar(1:length(valexp),100*valexp,100*storeEXPstd(iexp))
hold off
title('% of grasp outside box /user Experts')
xlabel('Surgeon #')
ylabel('% of grasp outside box')