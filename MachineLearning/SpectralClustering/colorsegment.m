function labels = colorsegment(I,k);
% 
% This function segments an image using color and proximity cues.
% 
% I is the input image
% k is the number of groups
% labels is the output partition image.

[N1,N2,N3]=size(I);
if N3 == 1;
  sigmac=0.3;
else
  sigmac=0.03;
end

sigmap=70;

n=25;
threshold=0.0;
niter=5;
%nsamp=min(max((k*0.001*N1*N2),50),100);
nsamp=100;
labels=simple(I,k,sigmac,sigmap,nsamp,n,threshold,niter);