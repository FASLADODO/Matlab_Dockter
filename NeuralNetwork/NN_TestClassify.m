function accuracy = NN_TestClassify(Data, NN, Classes)

classify = [];
for tt = 1:length(Data);
    %test classification
    outputNode = NN_Online(NN, Data{tt}.inputs);
    error = [];
    for ii = 1:length(Classes)
        error(ii) = norm(Classes(ii,:) - outputNode{end} );
    end
    [val,idx] = min(error);
    %store classes
    classify = [classify; Data{tt}.class, idx];
end

temp = [];
for ii = 1:length(classify)
    %Check over time if error2 was less than error 1
    temp = [temp; classify(ii,1) == classify(ii,2)];
end

accuracy = sum(temp)/length(temp); %percent of correct

end 