function labels = simple(I,k,sigmac,sigmap,nsamp,n,t,niter)
% function labels = ncut_segment(I,sigmap,sigmac)
%
% This function returns a partition image for the given image I
% It takes two parameters, 
%
%      k : Number of groups to extract
% sigmac : color change sensitivity
% sigmap : proximity sensitivity
%  nsamp : number of samples to be used
%      n : number of iterations of kmeans
%      t : k-means termination threshold
%  niter : number of randomized restarts for k-means

% how do you decide on nsamp

[N1,N2,N3]=size(I);

if max(I(:)) > 1.0;
  cdata=reshape(double(I)/255,[N1*N2 N3]);
else
  cdata=reshape(double(I),[N1*N2 N3]);
end

[X,Y]=meshgrid(1:N2,1:N1);
pdata=[X(:) Y(:)];
N=N1*N2;
%nsamp= max(round(k*0.01*N),50);

ind_all=randperm(N);
ind_samp=ind_all(1:nsamp);
ind_rest=ind_all(nsamp+1:end);


d2Ac=pairdist(cdata(ind_samp,:),cdata(ind_samp,:),'L2');
d2Ap=pairdist(pdata(ind_samp,:),pdata(ind_samp,:),'L2');
d2Bc=pairdist(cdata(ind_samp,:),cdata(ind_rest,:),'L2');
d2Bp=pairdist(pdata(ind_samp,:),pdata(ind_rest,:),'L2');

A=exp(-d2Ac.^2/(2*sigmac^2) -(d2Ap.^2/(2*sigmap^2))); % choose your sigma 
B=exp(-d2Bc.^2/(2*sigmac^2) -d2Bp.^2/(2*sigmap^2));

% the joys of nystrom
[V,ss]=nystrom_ncut(A,B);

% depermute the nsytrom result
V(ind_all,:)=V;
V1=V(:,1:k); % choose the k eigenvectors
Vnormalized=V1./repmat(sqrt(sum(V1.^2,2)),1,k);
[centers,labels,error]=kmwrapper(Vnormalized',k,n,t,niter);
labels=reshape(labels,[N1 N2]);

