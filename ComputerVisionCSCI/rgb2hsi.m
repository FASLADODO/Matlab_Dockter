function [hsiImage] = rgb2hsi(inimage)
if length(size(inimage))~=3
    error('rgb image needed for hsi conversion');
end
inimage=im2double(inimage);
rmat=inimage(:,:,1);
gmat=inimage(:,:,2);
bmat=inimage(:,:,3);
Hue=acos((0.5*((rmat-gmat)+(rmat-bmat)))./((sqrt((rmat-gmat).^2+(rmat-bmat).*(gmat-bmat)))+eps));
Hue(bmat>gmat)=2*pi-Hue(bmat>gmat);
Hue=Hue/(2*pi);
Sat=1-3.*(min(min(rmat,gmat),bmat))./(rmat+gmat+bmat+eps);
Int=(rmat+gmat+bmat)/3;
hsiImage=cat(3,Hue,Sat,Int);

end