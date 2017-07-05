%emboss
%apply emboss filter to image
%%

picnamestring = 'apple.jpg';
inputpic = imread(picnamestring);
imshow(inputpic)

rpic = inputpic(:,:,1);
gpic = inputpic(:,:,2);
bpic = inputpic(:,:,3);

mask = [-2,-1,0;-1,1,1;0,1,2];

rconv = filter2(mask,rpic);
gconv = filter2(mask,gpic);
bconv = filter2(mask,bpic);

newim = cat(3,rconv,gconv,bconv);
imshow(newim)

%%