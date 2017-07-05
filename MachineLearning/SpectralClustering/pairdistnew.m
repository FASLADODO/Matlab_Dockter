function d=pairdistnew(X,Y,dist_fun);
% d=pairdistnew(X,Y,dist_fun);
%
% finds the distance between all pairs of (row) vectors in arrays X and Y
%
% dist_fun can be 'L1', 'L2', 'Linf', 'chisq'
% faster vectorized versions 
%
% Sameer sagarwal@cs.ucsd.edu

[M1,N1]=size(X);
[M2,N2]=size(Y);

d=zeros(M1,M2);
if N1~=N2
   error('X and Y need to have the same number of columns')
end

% I think its more efficent to loop over the coordinates instead of vectorizing it
% since the number of coordinates is fairly small, compared to the original arrays
% this does not impose a significant loop overhead, plus it removes the need for the multiple
% copies of the array

if strcmp(dist_fun,'L1')
    for i = 1:N1
        T1=repmat(X(:,i),1,M2);
        T2=repmat(Y(:,i),1,M1)';
        d=d+abs(T1-T2); 
    end

elseif strcmp(dist_fun,'Linf')
    for i = 1:N1
        T1=repmat(X(:,i),1,M2);
        T2=repmat(Y(:,i),1,M1)';
 %       size(d);
 %       size(abs(T1-T2))
        d=max(d,abs(T1-T2)); 
    end

elseif strcmp(dist_fun,'chisq')
    for i = 1:N1
        T1=repmat(X(:,i),1,M2);
        T2=repmat(Y(:,i),1,M1)';
        d=d+(T2-T1).^2./(T1+T2+eps);
    end
    d=d.*0.5;
    
elseif strcmp(dist_fun,'L2')
    for i = 1:N1
        T1=repmat(X(:,i),1,M2);
        T2=repmat(Y(:,i),1,M1)';
        d=d+(T2-T1).^2;
    end
    d=sqrt(d);
elseif strcmp(dist_fun,'KS')
    disp('wow')
    for i = 1:M1
        x1=X(i,:);
        for j = 1:M2
            y1=Y(j,:);
            d(i,j)=kstest3(x1',y1');
        end
    end
            
            
        
end
