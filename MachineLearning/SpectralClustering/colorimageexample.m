I=imread('S_beer.jpeg');

figure(1);
imagesc(I);
title('Original Image');
axis image;
axis off;

[N1,N2,N3]=size(I);
k=5;
sigmac=0.3;
sigmap=70;
n=25;
threshold=0.0;
niter=5;
nsamp=100;

labels=simple(I,k,sigmac,sigmap,nsamp,n,threshold,niter);

figure(2);
imagesc(labels);
axis image;
axis off;
title('Partition Image');

figure(3);
imagesc(meanfill(I,Iseg));
axis image;
axis off;
title('Average Color Image');