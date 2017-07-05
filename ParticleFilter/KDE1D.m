function [ prob, x ] = KDE1D( dataset, kernalWidth )
%KDE_1D solves the KDE value for a 1d dataset
%Using Parzen windows
%http://research.cs.tamu.edu/prism/lectures/pr/pr_l7.pdf
%http://www.mathworks.com/help/stats/kernel-distribution.html
%f(x)=1/nh*sum(K(x?xi))

    lowerB = min(dataset) - range(dataset)*0.1;
    upperB = max(dataset) + range(dataset)*0.1;
    
    
    pd = fitdist(dataset,'Kernel','BandWidth',kernalWidth);
    x = lowerB:.1:upperB;
    y = pdf(pd,x);
    
    %plot
%     figure(1);
%     plot(x,y,'k-','LineWidth',2);
%     
%     figure(2);
%     hist(dataset);
    
    prob = y;
end

