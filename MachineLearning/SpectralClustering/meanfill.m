function value = meanfill(I,labels);
[N1,N2,N3]=size(I);
data=reshape(I,[N1*N2 N3]);
nmax=max(labels(:));
nmin=min(labels(:));
ldata=labels(:);
value=zeros(N1*N2,N3);
for i = nmin:nmax;
  slice=find(labels==i);
  M1=length(slice);
  value(slice,:)=repmat(round(mean(data(slice,:))),M1,1);
end

value=uint8(reshape(value,[N1 N2 N3]));