function MCC = confusionMatMCC(cm)

%cm = confusionmat(known, predicted);

TP = cm(1,1);
FN = cm(2,1);
FP = cm(1,2);
TN = cm(2,2);

MCC = (TP*TN - FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));

% return zero if divide by zero occurs (very crappy classifier)
if isequalwithequalnans(MCC, NaN)
    MCC = 0;
end