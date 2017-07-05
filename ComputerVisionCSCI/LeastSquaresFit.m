function [slope, intercept] = LeastSquaresFit(Xpoints, Ypoints)
%%Extra Credit 4
%%Rodney Dockter
%%reused for hough transform

if length(Xpoints) ~= length(Ypoints) 
    error('Coordinates should be same size')
end
X=zeros(length(Xpoints),2);
X(:,2)=Xpoints';
X(:,1)=ones(length(Xpoints),1);
Y=Ypoints';

B=(transpose(X)*X\transpose(X))*Y;
slope = B(2,1);
intercept = B(1,1);
end

