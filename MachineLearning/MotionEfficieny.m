function [D,X1,X2,params] = MotionEfficieny(Data)
    %distance to minimum path line at each point in Data
    %http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
    %Data should be an xyz matrix with rows as samples
    
    [NN,SS] = size(Data);
    %window = round(0.03*NN);
    window = 1;
    if(window < 1)
       window = 1; 
    end
    
    %get sorted data
    %DS = SortNorm(Data);
    
%     Data(1:window,:)
%     Data(NN-window+1:NN,:)
    %get point1;
    %X1 = mean(Data(1:window,:));
    X1 = Data(1:window,:);
    
    %get point2
    %X2 = mean(Data(NN-window+1:NN,:));
    X2 = Data(NN-window+1:NN,:);
    
    %Equation 10 numerator
    shiftx1 = bsxfun(@minus,Data,X1);
    shiftx2 = bsxfun(@minus,Data,X2);
    numerator = NormRowWise( cross(shiftx1,shiftx2) );
    
    %Equation 10 denominator (scale by the overall length of the line)
    denominator = NormRowWise(X2-X1);
    
    %result
    D = numerator ./ denominator;
    
    %also get line params for kicks
    params = [];
    if(false)
        xx = [X1; X2];
        params = pinv([xx(:,[1,2]),1;1])*xx(:,[3]);
    end
end