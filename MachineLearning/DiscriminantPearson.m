function [DPWEIGHT,RANK] = DiscriminantPearson(X,Labels)
%use pearson correlation idea to determine between class seperability
%returns DP = seperability for each state

%siget sizes and number of labels
[NN,SS] = size(X);
cslist = unique(Labels);

%scale by range of data
XMIN = min(X);
XMAX = max(X);
XDIFF = XMAX - XMIN;
X = bsxfun(@rdivide,X,XDIFF);

%loop through all classes
for ii = 1:length(cslist);
    DON = X(Labels == cslist(ii),:);
    DOFF = X(Labels ~= cslist(ii),:);
    for jj = 1:SS
       %get state data
        xon = DON(:,jj); 
        xoff = DOFF(:,jj); 
        xall = X(:,jj); %both classes
        
        %zero mean shift
        Eon = xon - mean(xon);
        Eoff = xoff - mean(xoff);
        Eall = xall - mean(xall); %both classes

        %compute dis-correlation?
        denom = ( sqrt(sum(Eon.^2))*sqrt(sum(Eoff.^2)) );
        if(denom > eps)
            DP(ii,jj) =  (sqrt(sum(Eall.^2))) / denom;
        else
            
        end
    end
end

DPWEIGHT = mean(DP); %just give mean values for each state

%sort em
[~,RANK] = sort(DPWEIGHT,'descend');

end