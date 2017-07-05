function e = ClassificationError(Class_Estimate,Class_Actual)
    
    err = Class_Estimate ~= Class_Actual;
    
    e = sum(err)/length(err);
end