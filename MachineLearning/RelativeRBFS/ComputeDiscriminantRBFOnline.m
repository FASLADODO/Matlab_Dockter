function [PROnline] = ComputeDiscriminantRBFOnline(Data_Online,Data_All,bandwidth)
%Data_Online is matrix of current data (rows are samples)
%Data_All is struct for training data with each entry for different class

    clist = 1:length(Data_All);
    
    for cc = clist
        nn = size(Data_All{cc},1);
        %get pdists
        dist = pdist2(Data_Online,Data_All{cc});
        distinternal = pdist2(Data_All{cc},Data_All{cc});
        %get rbfs
        %rw = exp( -(bandwidth*dist).^2 );
        rw = exp( -(dist.^2)./(2*bandwidth^2) );
        bi = exp( -(distinternal.^2)./(2*bandwidth^2) );
        %get scale within
        internal_us = sum(bi,2);
        scalei = ScaleRBF(Data_All{cc},internal_us);
        %get within and between class
        pon{cc} = sum(rw,2)/scalei;
    end
    PROnline{1} = ComputeRBFDifference(pon{1},pon{2});
    PROnline{2} = ComputeRBFDifference(pon{2},pon{1});
end