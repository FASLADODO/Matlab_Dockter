function windowData = DataInWindow(Data,CPoint,windowSize)
%Return subset of Data that is within windowSize of Cpoint
%Data: matrix with columns as states, rows as samples
%CPoint: 1xn vector indicating the center point of the window (n=#states)
%windowSize: scalar value which indicates the size of windows

    %start function
    [NN,SS] = size(Data);
    
    if(length(CPoint) ~= SS)
       error('Center Point size does not match # states in Data') 
    end

    orgn = repmat(CPoint,NN,1);

    shiftedData = Data - orgn;

    normShift = sqrt(sum(abs(shiftedData.^2),2));

    idx = find(normShift < windowSize);

    windowData = Data(idx,:);

    %Se Fin
end